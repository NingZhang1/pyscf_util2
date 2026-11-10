#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""MRPT2 积分导出：优化版 vs 旧版 的正确性与耗时对比。

用法::

    python examples/16-MRPT2-opt.py              # 默认 cc-pVTZ（约 20-30 秒）
    python examples/16-MRPT2-opt.py cc-pVQZ      # 115 个轨道，差异更明显（约 2-4 分钟）
    python examples/16-MRPT2-opt.py cc-pVDZ      # 快速冒烟测试（约 5 秒）

本脚本做两件事：

1. **正确性**：把旧实现（integral_MRPT2.py / integral_MRPT2_incore_fast.py）
   与新实现（integral_MRPT2_opt.py / integral_MRPT2_binary.py）在同一分子上
   逐项对比并打印 PASS/FAIL：

   * incore / outcore  -> 输出文件 **逐字节相同**（md5 相等）
   * fast              -> 8-fold 键集合相同、逐键积分值差 <= 1e-10、无重复记录
   * binary            -> 轨道指标/积分值与 ASCII 一致；outcore 与 incore 的
                          二进制文件逐字节相同；转回 ASCII 后能被
                          pyscf.tools.fcidump.read 正确读回，并与 pyscf 自己的
                          fcidump.from_mo 参考结果一致（H1/H2/ECORE）

2. **性能**：对每一条路径分别计时并打印加速比与文件体积。
   所有计时都包含 AO -> MO 变换（这是新旧实现共有的部分）。

若任何一项检查失败，脚本以退出码 1 结束。
"""

import contextlib
import hashlib
import io
import os
import re
import shutil
import sys
import tempfile
import time
from collections import defaultdict

import numpy as np
import pyscf
from pyscf import gto, scf, tools

from pyscf_util.Integrals import integral_MRPT2 as old_incore_mod
from pyscf_util.Integrals import integral_MRPT2_incore_fast as old_fast_mod
from pyscf_util.Integrals import integral_MRPT2_opt as opt_mod
from pyscf_util.Integrals import integral_MRPT2_binary as bin_mod


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def key8(a, b, c, d):
    """8-fold 对称键（pair 形式，与 FCIDUMP 的 (pq|rs) 约定一致）。"""
    ij = (a, b) if a >= b else (b, a)
    kl = (c, d) if c >= d else (d, c)
    return (ij, kl) if ij >= kl else (kl, ij)


def parse_fcidump(path):
    """解析 FCIDUMP：返回 (records, h1e, ecore, orbsym, norb, nelec)。

    ``records`` 是 ``{8-fold 键: [值, ...]}``；重复记录会保留下来，便于检测
    旧 fast 版本中的冗余积分。
    """
    records = defaultdict(list)
    h1e = {}
    ecore = None
    orbsym = None
    norb = None
    nelec = None
    with open(path) as fin:
        lines = fin.read().splitlines()
    end = 0
    for end, line in enumerate(lines):
        if "&END" in line:
            break
    head = "\n".join(lines[: end + 1])
    match = re.search(r"NORB\s*=\s*(\d+)", head)
    if match:
        norb = int(match.group(1))
    match = re.search(r"NELEC\s*=\s*(\d+)", head)
    if match:
        nelec = int(match.group(1))
    match = re.search(r"ORBSYM\s*=\s*([0-9,\s]+)", head)
    if match:
        orbsym = [int(x) for x in re.split(r"[,\s]+", match.group(1).strip()) if x]
    for line in lines[end + 1 :]:
        fields = line.split()
        if len(fields) < 5:
            continue
        value = float(fields[0])
        i, j, k, l = (int(x) for x in fields[1:5])
        if k == 0:
            if j == 0:
                ecore = value
            else:
                h1e[(i, j)] = value
            continue
        records[key8(i, j, k, l)].append(value)
    return records, h1e, ecore, orbsym, norb, nelec


def md5_of(path):
    with open(path, "rb") as fin:
        return hashlib.md5(fin.read()).hexdigest()


def timed(call, *args, **kwargs):
    """执行 ``call(*args)`` 并返回 ``(result, elapsed_seconds)``。"""
    start = time.time()
    result = call(*args, **kwargs)
    return result, time.time() - start


def quiet(call, *args, **kwargs):
    """执行调用但吞掉旧模块打印的进度信息（保持计时输出整齐）。"""
    with contextlib.redirect_stdout(io.StringIO()):
        return timed(call, *args, **kwargs)


def paired_best(repeat, old_call, old_args, new_call, new_args):
    """交替计时旧/新调用 ``repeat`` 次，各取最小值，返回 (t_old, t_new, result)。

    交替 + 取最小值可以显著降低共享机器上的计时抖动（单次计时常见 ±15%）。
    """
    best_old = float("inf")
    best_new = float("inf")
    result = None
    for _ in range(max(1, repeat)):
        _, elapsed = quiet(old_call, *old_args)
        best_old = min(best_old, elapsed)
        result, elapsed = timed(new_call, *new_args)
        best_new = min(best_new, elapsed)
    return best_old, best_new, result


class Report(object):
    """收集耗时行与 PASS/FAIL 检查。"""

    def __init__(self):
        self.rows = []
        self.checks = []

    def add(self, task, t_old, t_new, size_old, size_new, note=""):
        self.rows.append((task, t_old, t_new, size_old, size_new, note))

    def check(self, label, ok, detail=""):
        self.checks.append((label, bool(ok), detail))
        print("  [%s] %-48s %s" % ("OK  " if ok else "FAIL", label, detail))

    def print_timings(self):
        print("")
        print("=" * 102)
        print(
            "%-36s %10s %10s %9s %16s" % ("任务", "旧版", "新版", "加速比", "文件体积")
        )
        print("-" * 102)
        for task, t_old, t_new, size_old, size_new, note in self.rows:
            speedup = t_old / t_new if t_new else float("inf")
            if size_old:
                size = "%.1fMB -> %.1fMB" % (size_old / 2.0**20, size_new / 2.0**20)
            else:
                size = "%.1fMB" % (size_new / 2.0**20)
            print(
                "%-36s %9.2fs %9.2fs %8.2fx %16s" % (task, t_old, t_new, speedup, size)
            )
            if note:
                print("        " + note)
        print("-" * 102)

    def ok(self):
        return all(flag for _, flag, _ in self.checks)


def _binary_records(path):
    """返回二进制文件的 (1-based 8-fold 键集合, {键: 值})。"""
    data = bin_mod.load_mrpt2_binary(path)
    keys = {}
    for record in data["eri"]:
        k = key8(
            int(record["p"]) + 1,
            int(record["q"]) + 1,
            int(record["r"]) + 1,
            int(record["s"]) + 1,
        )
        keys[k] = float(record["v"])
    return set(keys), keys, data


# ---------------------------------------------------------------------------
# main demo
# ---------------------------------------------------------------------------


def main(basis):
    mol = gto.M(
        verbose=0,
        atom="O 0 0 0; H 0 0.757 0.586; H 0 -0.757 0.586",
        basis=basis,
        symmetry="C2v",
        spin=0,
        charge=0,
    )
    mol.build()
    mf = scf.RHF(mol).run()
    mo = mf.mo_coeff
    nfzc, nact = 1, 8
    nvir = mol.nao - nfzc - nact
    nmo = mol.nao
    tol = 1e-10

    print("=" * 102)
    print("MRPT2 积分导出：旧版 vs 优化版")
    print("  分子      : H2O / %s (C2v)" % basis)
    print(
        "  轨道      : nao = nmo = %d, nfzc = %d, nact = %d, nvir = %d"
        % (nmo, nfzc, nact, nvir)
    )
    print("  阈值      : tol = %g" % tol)
    print("  计时说明  : 所有时间都包含 AO -> MO 积分变换（新旧实现共有）")
    print("=" * 102)

    tmpdir = tempfile.mkdtemp(prefix="mrpt2_opt_example_")
    report = Report()
    try:
        # ------------------------------------------------------------------
        # 1. incore: 逐字节相同
        # ------------------------------------------------------------------
        print("")
        print("[1/4] incore 版本（fcidump_mrpt2 -> fcidump_mrpt2_opt）")
        f_old_incore = os.path.join(tmpdir, "incore_old.txt")
        f_new_incore = os.path.join(tmpdir, "incore_new.txt")
        t_old_incore, t_new_incore, stats_incore = paired_best(
            2,
            old_incore_mod.fcidump_mrpt2,
            (mol, mf, mo, nfzc, nact, nvir, f_old_incore, tol),
            opt_mod.fcidump_mrpt2_opt,
            (mol, mf, mo, nfzc, nact, nvir, f_new_incore, tol),
        )
        incore_same = md5_of(f_old_incore) == md5_of(f_new_incore)
        report.check(
            "incore 输出逐字节相同",
            incore_same,
            "md5 %s, %d 条记录" % (md5_of(f_new_incore)[:8], stats_incore["nrecords"]),
        )
        rec_incore, h1_incore, e_incore, _, _, _ = parse_fcidump(f_old_incore)
        report.add(
            "incore (ASCII)",
            t_old_incore,
            t_new_incore,
            os.path.getsize(f_old_incore),
            os.path.getsize(f_new_incore),
        )

        # ------------------------------------------------------------------
        # 2. outcore: 逐字节相同 + h5 复用
        # ------------------------------------------------------------------
        print("")
        print(
            "[2/4] outcore 版本（fcidump_mrpt2_outcore -> fcidump_mrpt2_outcore_opt）"
        )
        f_old_outcore = os.path.join(tmpdir, "outcore_old.txt")
        f_new_outcore = os.path.join(tmpdir, "outcore_new.txt")
        _, t_old_outcore = quiet(
            old_incore_mod.fcidump_mrpt2_outcore,
            mol,
            mf,
            mo,
            nfzc,
            nact,
            nvir,
            f_old_outcore,
            tol,
        )
        stats_outcore, t_new_outcore = timed(
            opt_mod.fcidump_mrpt2_outcore_opt,
            mol,
            mf,
            mo,
            nfzc,
            nact,
            nvir,
            f_new_outcore,
            tol,
        )
        outcore_same = md5_of(f_old_outcore) == md5_of(f_new_outcore)
        report.check(
            "outcore 输出逐字节相同",
            outcore_same,
            "md5 %s, %d 条记录"
            % (md5_of(f_new_outcore)[:8], stats_outcore["nrecords"]),
        )
        report.add(
            "outcore (ASCII)",
            t_old_outcore,
            t_new_outcore,
            os.path.getsize(f_old_outcore),
            os.path.getsize(f_new_outcore),
        )

        h5file = os.path.join(tmpdir, "eri_mo.h5")
        f_reuse_a = os.path.join(tmpdir, "reuse_a.txt")
        f_reuse_b = os.path.join(tmpdir, "reuse_b.txt")
        opt_mod.fcidump_mrpt2_outcore_opt(
            mol, mf, mo, nfzc, nact, nvir, f_reuse_a, tol, erifile=h5file
        )
        _, t_reuse = timed(
            opt_mod.fcidump_mrpt2_outcore_opt,
            mol,
            mf,
            mo,
            nfzc,
            nact,
            nvir,
            f_reuse_b,
            tol,
            erifile=h5file,
        )
        report.check(
            "outcore 复用已有 MO 积分文件（不重做 AO->MO）",
            md5_of(f_reuse_a) == md5_of(f_reuse_b),
            "第二次调用只需 %.2fs" % t_reuse,
        )
        report.add(
            "outcore 复用 MO 积分 (ASCII)",
            t_new_outcore,
            t_reuse,
            os.path.getsize(f_new_outcore),
            os.path.getsize(f_reuse_b),
            note="上：首次含 AO->MO 变换；下：复用 h5，只剩记录枚举 + 写文件",
        )

        # ------------------------------------------------------------------
        # 3. fast: 键集合相同、值一致、去重
        # ------------------------------------------------------------------
        print("")
        print(
            "[3/4] fast 版本（fcidump_mrpt2_incore_fast -> fcidump_mrpt2_incore_fast_opt）"
        )
        f_old_fast = os.path.join(tmpdir, "fast_old.txt")
        f_new_fast = os.path.join(tmpdir, "fast_new.txt")
        _, t_old_fast = quiet(
            old_fast_mod.fcidump_mrpt2_incore_fast,
            mol,
            mf,
            mo,
            nfzc,
            nact,
            nvir,
            f_old_fast,
            tol,
        )
        stats_fast, t_new_fast = timed(
            opt_mod.fcidump_mrpt2_incore_fast_opt,
            mol,
            mf,
            mo,
            nfzc,
            nact,
            nvir,
            f_new_fast,
            tol,
        )
        rec_fast_old, h1_fast_old, e_fast_old, sym_fast_old, _, _ = parse_fcidump(
            f_old_fast
        )
        rec_fast_new, h1_fast_new, e_fast_new, sym_fast_new, _, _ = parse_fcidump(
            f_new_fast
        )
        n_lines_old = sum(len(v) for v in rec_fast_old.values())
        n_lines_new = sum(len(v) for v in rec_fast_new.values())
        max_delta = max(
            (
                abs(rec_fast_old[k][0] - rec_fast_new[k][0])
                for k in set(rec_fast_old) & set(rec_fast_new)
            ),
            default=0.0,
        )
        report.check(
            "fast 8-fold 键集合相同",
            set(rec_fast_old) == set(rec_fast_new),
            "%d 个键；旧版写 %d 条（含重复），新版写 %d 条"
            % (len(rec_fast_old), n_lines_old, n_lines_new),
        )
        report.check(
            "fast 逐键积分值一致（<=1e-10）",
            max_delta <= 1e-10,
            "max|Δ| = %.2e" % max_delta,
        )
        report.check(
            "fast 新版无重复记录",
            n_lines_new == len(rec_fast_new),
            "重复记录 %d 条（旧版 %d 条）"
            % (n_lines_new - len(rec_fast_new), n_lines_old - len(rec_fast_old)),
        )
        report.check(
            "fast h1e / ECORE / ORBSYM 一致",
            h1_fast_old == h1_fast_new
            and abs(e_fast_old - e_fast_new) < 1e-12
            and sym_fast_old == sym_fast_new,
            "h1e %d 项, ECORE 差 %.1e"
            % (len(h1_fast_new), abs(e_fast_old - e_fast_new)),
        )
        report.add(
            "fast (ASCII)",
            t_old_fast,
            t_new_fast,
            os.path.getsize(f_old_fast),
            os.path.getsize(f_new_fast),
        )

        # ------------------------------------------------------------------
        # 4. binary: 往返一致 / outcore==incore / 与 pyscf 参考一致
        # ------------------------------------------------------------------
        print("")
        print("[4/4] 二进制格式（integral_MRPT2_binary）")
        f_bin_incore = os.path.join(tmpdir, "incore.bin")
        f_bin_outcore = os.path.join(tmpdir, "outcore.bin")
        stats_bin, t_bin_incore = timed(
            bin_mod.dump_mrpt2_binary,
            mol,
            mf,
            mo,
            nfzc,
            nact,
            nvir,
            f_bin_incore,
            tol,
            "mrpt2",
        )
        _, t_bin_outcore = timed(
            bin_mod.dump_mrpt2_binary_outcore,
            mol,
            mf,
            mo,
            nfzc,
            nact,
            nvir,
            f_bin_outcore,
            tol,
            "mrpt2",
        )
        report.check(
            "binary outcore == incore（逐字节）",
            md5_of(f_bin_incore) == md5_of(f_bin_outcore),
            "%d 条记录, %.1f MB" % (stats_bin["nint"], stats_bin["bytes"] / 2.0**20),
        )

        bin_keys, bin_values, data = _binary_records(f_bin_incore)
        header = bin_mod.read_mrpt2_binary_header(f_bin_incore)
        report.check(
            "binary 记录集合与 incore ASCII 一致",
            bin_keys == set(rec_incore),
            "%d 个 8-fold 键, NINT=%d" % (len(bin_keys), header["nint"]),
        )
        bin_delta = max(
            (abs(bin_values[k] - rec_incore[k][0]) for k in bin_keys), default=0.0
        )
        report.check(
            "binary 积分值与 incore ASCII 一致",
            bin_delta <= 1e-12,
            "max|Δ| = %.2e（ASCII 只保留 16 位有效数字）" % bin_delta,
        )
        h1_delta = max(
            abs(data["h1e"][i - 1, j - 1] - value)
            for (i, j), value in h1_incore.items()
        )
        report.check(
            "binary h1e / ECORE 一致",
            h1_delta <= 1e-12 and abs(data["nuc"] - e_incore) <= 1e-12,
            "max|Δh1e| = %.2e, |ΔECORE| = %.2e"
            % (h1_delta, abs(data["nuc"] - e_incore)),
        )
        report.check(
            "binary 体积小于 ASCII 的 55%",
            stats_bin["bytes"] < 0.55 * os.path.getsize(f_old_incore),
            "%.1f MB vs %.1f MB（16 字节/条：4×uint16 + float64）"
            % (stats_bin["bytes"] / 2.0**20, os.path.getsize(f_old_incore) / 2.0**20),
        )
        report.add(
            "incore (binary, 同一套积分)",
            t_old_incore,
            t_bin_incore,
            os.path.getsize(f_old_incore),
            stats_bin["bytes"],
            note="上：旧版写 ASCII；下：新版写 binary（都包含 AO->MO 变换）",
        )

        f_bin_fast = os.path.join(tmpdir, "fast.bin")
        stats_bin_fast, t_bin_fast = timed(
            bin_mod.dump_mrpt2_binary_fast,
            mol,
            mf,
            mo,
            nfzc,
            nact,
            nvir,
            f_bin_fast,
            tol,
        )
        fast_bin_keys, _, _ = _binary_records(f_bin_fast)
        report.check(
            "binary fast 记录集合与 fast ASCII 一致",
            fast_bin_keys == set(rec_fast_old),
            "%d 条, %.1f MB"
            % (stats_bin_fast["nint"], stats_bin_fast["bytes"] / 2.0**20),
        )
        report.add(
            "fast (binary)",
            t_old_fast,
            t_bin_fast,
            os.path.getsize(f_old_fast),
            stats_bin_fast["bytes"],
            note="上：旧版写 ASCII；下：新版写 binary（都包含 AO->MO 变换）",
        )

        # full 模式与 pyscf 自己的 FCIDUMP 参考对比
        if nmo <= 80:
            f_bin_full = os.path.join(tmpdir, "full.bin")
            stats_full, t_full = timed(
                bin_mod.dump_mrpt2_binary,
                mol,
                mf,
                mo,
                nfzc,
                nact,
                nvir,
                f_bin_full,
                tol,
                "full",
            )
            f_ref = os.path.join(tmpdir, "pyscf_reference.txt")
            _, t_ref = timed(tools.fcidump.from_mo, mol, f_ref, mo, tol=tol)
            f_conv = os.path.join(tmpdir, "full_converted.txt")
            _, t_conv = timed(bin_mod.mrpt2_binary_to_fcidump, f_bin_full, f_conv)
            ref = tools.fcidump.read(f_ref, molpro_orbsym=False)
            got = tools.fcidump.read(f_conv, molpro_orbsym=False)
            h1_full = np.abs(ref["H1"] - got["H1"]).max()
            h2_full = np.abs(ref["H2"] - got["H2"]).max()
            e_full = abs(ref["ECORE"] - got["ECORE"])
            report.check(
                "binary full 与 pyscf fcidump.from_mo 一致",
                h1_full <= 1e-12 and h2_full <= 1e-12 and e_full <= 1e-12,
                "H1 %.1e, H2 %.1e, ECORE %.1e（%d 条记录）"
                % (h1_full, h2_full, e_full, stats_full["nint"]),
            )
            report.check(
                "binary -> ASCII 可被 pyscf.tools.fcidump.read 读取",
                ref["NORB"] == got["NORB"] and ref["NELEC"] == got["NELEC"],
                "NORB=%d, NELEC=%d（转换耗时 %.2fs）"
                % (got["NORB"], got["NELEC"], t_conv),
            )
            report.add(
                "full 积分：pyscf ASCII vs binary",
                t_ref,
                t_full,
                os.path.getsize(f_ref),
                stats_full["bytes"],
                note="上：pyscf 自己写 ASCII；下：二进制写出（同样的积分集合）",
            )
        else:
            print(
                "  [SKIP] %-48s %s"
                % (
                    "binary full vs pyscf fcidump.from_mo",
                    "nao=%d 时全空间积分为 nmo^4/24 ~ %d 条，"
                    "pyscf 的 ASCII 参考过慢；该检查由小体系测试覆盖"
                    % (nmo, int(round(nmo**4 / 24.0))),
                )
            )

        # 读取速度
        _, t_load = timed(bin_mod.load_mrpt2_binary, f_bin_incore)
        ascii_read_start = time.time()
        tools.fcidump.read(f_old_incore, molpro_orbsym=False)
        t_ascii_read = time.time() - ascii_read_start
        report.check(
            "binary 读取快于 ASCII 解析",
            t_load < t_ascii_read,
            "load %.3fs vs fcidump.read %.3fs" % (t_load, t_ascii_read),
        )
        report.add(
            "读取：ASCII 解析 vs binary load",
            t_ascii_read,
            t_load,
            os.path.getsize(f_old_incore),
            stats_bin["bytes"],
        )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    report.print_timings()
    ascii_rows = [
        row
        for row in report.rows
        if row[0] in ("incore (ASCII)", "outcore (ASCII)", "fast (ASCII)")
    ]
    if len(ascii_rows) == 3:
        sum_old = sum(row[1] for row in ascii_rows)
        sum_new = sum(row[2] for row in ascii_rows)
        print(
            "三条 ASCII 路径合计: 旧版 %.2fs -> 新版 %.2fs (%.2fx)"
            % (sum_old, sum_new, sum_old / sum_new if sum_new else float("inf"))
        )
        print(
            "（incore/outcore 为交替两次取最小值；其余为单次计时，"
            "共享机器上单次值可能有 ±10-15% 抖动）"
        )
    print("")
    if report.ok():
        print("全部检查通过：优化版输出与旧版一致，二进制格式可正确往返。")
        return 0
    print("存在失败项，请检查上面的 [FAIL] 输出。")
    return 1


if __name__ == "__main__":
    _basis = sys.argv[1] if len(sys.argv) > 1 else "cc-pVTZ"
    sys.exit(main(_basis))
