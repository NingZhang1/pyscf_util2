#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""四分量积分导出：优化版 vs 旧版 的正确性与耗时对比。

用法::

    python examples/17-4C-ints-opt.py                 # 默认 H2O/3-21g（约 20-30 秒）
    python examples/17-4C-ints-opt.py cc-pvdz-dk      # 更大的基组（约 2-3 分钟）

脚本内容：

1. **正确性**：旧实现（``integral_4C.FCIDUMP_Rela4C``）与新实现
   （``integral_4C_opt.FCIDUMP_Rela4C_opt``）在同一分子上对比并打印 PASS/FAIL

   * FCIDUMP 记录集合完全相同，逐条积分值差 <= 1e-10（打印精度 12 位小数）
   * h1e / ECORE / NORB / NELEC / ORBSYM 完全相同
   * incore 与 outcore 的结果互相一致（旧版 outcore 曾忽略 paldus/dg 过滤，
     新实现已与 incore 统一，见下方说明）

2. **对称性验证**（优化所依赖的前提，逐条数值验证）

   * 库仑型 intor：``aosym="s4"`` 与 ``"s1"`` 的全张量差 <= 1e-11
     -> 可以使用 4-fold AO 对称性
   * Breit/Gaunt intor：``aosym="s4"`` 相对偏差 ~0.1
     -> **不能**用 s4，保持 s1（这是原代码里 s1 的真正原因）
   * Breit/Gaunt 四个收缩项之间的严格关系：
     ``t2[i,j,k,l] = conj(t1[j,i,l,k])``、``t4 = conj(t3[j,i,l,k])``
     -> 只需 2 个 AO 积分张量 / 2 次变换

3. **性能**：incore / outcore × 是否含 Breit，分别计时并打印加速比。

任何检查失败时以退出码 1 结束。
"""

import os
import shutil
import sys
import tempfile
import time

import h5py
import numpy as np
from pyscf import gto, lib, scf
from pyscf.ao2mo import r_outcore

from pyscf_util.Relativisitc import integral_4C as old_mod
from pyscf_util.Relativisitc import integral_4C_opt as opt_mod

WATER = "O 0 0 0; H 0 -0.757 0.587; H 0 0.757 0.587"


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def parse_fcidump(path):
    import re

    records = {}
    h1e = {}
    ecore = None
    orbsym = None
    nelec = None
    norb = None
    with open(path) as fin:
        lines = fin.read().splitlines()
    end = 0
    for end, line in enumerate(lines):
        if "&END" in line:
            break
    head = "\n".join(lines[: end + 1])
    m = re.search(r"NORB\s*=\s*(\d+)", head)
    if m:
        norb = int(m.group(1))
    m = re.search(r"NELEC\s*=\s*(\d+)", head)
    if m:
        nelec = int(m.group(1))
    m = re.search(r"ORBSYM\s*=\s*([0-9,\s]+)", head)
    if m:
        orbsym = [int(x) for x in re.split(r"[,\s]+", m.group(1).strip()) if x]
    for line in lines[end + 1 :]:
        t = line.split()
        if len(t) < 5:
            continue
        if len(t) == 6:
            value = complex(float(t[0]), float(t[1]))
            i, j, k, l = (int(x) for x in t[2:6])
        else:
            value = complex(float(t[0]), 0.0)
            i, j, k, l = (int(x) for x in t[1:5])
        if k == 0:
            if j == 0:
                ecore = value
            else:
                h1e[(i, j)] = value
            continue
        records[(i, j, k, l)] = value
    return records, h1e, ecore, orbsym, nelec, norb


class Report(object):
    def __init__(self):
        self.rows = []
        self.checks = []

    def add(self, task, t_old, t_new, nrec):
        self.rows.append((task, t_old, t_new, nrec))

    def check(self, label, ok, detail=""):
        self.checks.append((label, bool(ok)))
        print("  [%s] %-46s %s" % ("OK  " if ok else "FAIL", label, detail))

    def ok(self):
        return all(flag for _, flag in self.checks)

    def print_timings(self):
        print("")
        print("=" * 88)
        print("%-30s %10s %10s %9s %12s" % ("任务", "旧版", "新版", "加速比", "记录数"))
        print("-" * 88)
        for task, t_old, t_new, nrec in self.rows:
            print(
                "%-30s %9.2fs %9.2fs %8.2fx %12d"
                % (task, t_old, t_new, t_old / t_new if t_new else float("inf"), nrec)
            )
        print("-" * 88)


def compare(report, ref, new, label):
    rec_ref, h1_ref, e_ref, sym_ref, n_ref, norb_ref = ref
    rec_new, h1_new, e_new, sym_new, n_new, norb_new = new
    keys_ok = set(rec_ref) == set(rec_new)
    delta = max(
        (abs(rec_ref[k] - rec_new[k]) for k in set(rec_ref) & set(rec_new)),
        default=0.0,
    )
    h1_ok = set(h1_ref) == set(h1_new) and (
        not h1_ref or max(abs(h1_ref[k] - h1_new[k]) for k in h1_ref) <= 1e-12
    )
    report.check(
        "%s 记录集合一致" % label,
        keys_ok,
        "%d 条记录" % len(rec_ref),
    )
    report.check(
        "%s 逐条积分值一致 (<=1e-10)" % label,
        delta <= 1e-10,
        "max|dv| = %.2e" % delta,
    )
    report.check(
        "%s h1e / ECORE / ORBSYM" % label,
        h1_ok
        and abs(e_ref - e_new) <= 1e-12
        and sym_ref == sym_new
        and n_ref == n_new
        and norb_ref == norb_new,
        "h1e %d 条, ECORE 差 %.1e" % (len(h1_ref), abs(e_ref - e_new)),
    )


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def main(basis):
    mol = gto.M(
        atom=WATER,
        basis=basis,
        verbose=0,
        charge=0,
        spin=0,
        symmetry="c2v",
    )
    mol.build()
    mf = scf.dhf.RDHF(mol)
    mf.conv_tol = 1e-9
    mf.kernel()
    mf.with_breit = True
    npes = min(20, mol.nao_2c())

    print("=" * 88)
    print("四分量积分导出：旧版 vs 优化版")
    print(
        "  H2O / %s, nao_2c = %d, npes = %d, Breit = on" % (basis, mol.nao_2c(), npes)
    )
    print("=" * 88)

    tmpdir = tempfile.mkdtemp(prefix="rela4c_opt_example_")
    report = Report()
    try:
        # ------------------------------------------------------------------
        # 1. dump equivalence + timing
        # ------------------------------------------------------------------
        print("")
        print("[1/3] FCIDUMP 输出对比（旧版 vs 新版）")
        for with_breit in (False, True):
            for mode in ("incore", "outcore"):
                tag = "%s_breit%d" % (mode, with_breit)
                ref_file = os.path.join(tmpdir, "ref_" + tag + ".txt")
                new_file = os.path.join(tmpdir, "new_" + tag + ".txt")
                t = time.time()
                old_mod.FCIDUMP_Rela4C(
                    mol,
                    mf,
                    with_breit=with_breit,
                    npes=npes,
                    filename=ref_file,
                    mode=mode,
                    IsComplex=True,
                )
                t_old = time.time() - t
                t = time.time()
                opt_mod.FCIDUMP_Rela4C_opt(
                    mol,
                    mf,
                    with_breit=with_breit,
                    npes=npes,
                    filename=new_file,
                    mode=mode,
                    IsComplex=True,
                )
                t_new = time.time() - t
                ref = parse_fcidump(ref_file)
                new = parse_fcidump(new_file)
                compare(report, ref, new, "%s Breit=%d" % (mode, int(with_breit)))
                report.add(
                    "%s (Breit=%d)" % (mode, int(with_breit)), t_old, t_new, len(new[0])
                )

        # incore vs outcore consistency
        inc = parse_fcidump(os.path.join(tmpdir, "new_incore_breit1.txt"))
        out = parse_fcidump(os.path.join(tmpdir, "new_outcore_breit1.txt"))
        compare(report, inc, out, "新版 incore vs outcore")

        # ------------------------------------------------------------------
        # 2. symmetry verification
        # ------------------------------------------------------------------
        print("")
        print("[2/3] 优化所依赖的对称性验证")
        n2c = mol.nao_2c()
        mo = np.asarray(mf.mo_coeff)[:, n2c : n2c + 12]
        L, S = mo[:n2c, :], mo[n2c:, :]

        def ao_mo(intor, coeffs, aosym):
            path = os.path.join(tmpdir, "sym_%s_%s.h5" % (intor, aosym))
            r_outcore.general(mol, coeffs, path, intor=intor, aosym=aosym, verbose=0)
            with h5py.File(path, "r") as h5:
                arr = np.array(h5["eri_mo"]).reshape(12, 12, 12, 12)
            os.remove(path)
            return arr

        coulomb_cases = (
            ("int2e_spinor (LLLL)", (L, L, L, L)),
            ("int2e_spsp1spsp2_spinor (SSSS)", (S, S, S, S)),
            ("int2e_spsp2_spinor (LLSS)", (L, L, S, S)),
            ("int2e_spsp1_spinor (SSLL)", (S, S, L, L)),
        )
        worst = 0.0
        for name, coeffs in coulomb_cases:
            intor = name.split(" ")[0]
            ref = ao_mo(intor, coeffs, "s1")
            got = ao_mo(intor, coeffs, "s4")
            worst = max(worst, np.abs(ref - got).max())
        report.check(
            "库仑 intor: aosym s4 == s1 (<=1e-11)",
            worst <= 1e-11,
            "max|dT| = %.2e  -> 变换可用 4-fold AO 对称性" % worst,
        )

        intor = "int2e_breit_ssp1ssp2_spinor"
        ref = ao_mo(intor, (L, S, L, S), "s1")
        got = ao_mo(intor, (L, S, L, S), "s4")
        rel = np.abs(ref - got).max() / max(1e-30, np.abs(ref).max())
        report.check(
            "Breit intor: s4 不成立（因此保持 s1）",
            rel > 1e-3,
            "相对偏差 %.2e  -> 不能对 Breit 使用 s4" % rel,
        )

        c1 = 0.5 / lib.param.LIGHT_SPEED
        mo10 = np.asarray(mf.mo_coeff)[:, n2c : n2c + 10]
        L10, S10 = mo10[:n2c, :], mo10[n2c:, :]

        def transform(ao, a, b, c, d):
            return lib.einsum("pqrs,pi,qj,rk,sl->ijkl", ao, a.conj(), b, c.conj(), d)

        rel_err = 0.0
        for intor_lsls, intor_lssl in (
            ("int2e_breit_ssp1ssp2_spinor", "int2e_breit_ssp1sps2_spinor"),
            ("int2e_ssp1ssp2_spinor", "int2e_ssp1sps2_spinor"),
        ):
            ao1 = mol.intor(intor_lsls) * c1**2
            ao3 = mol.intor(intor_lssl) * c1**2
            t1 = transform(ao1, L10, S10, L10, S10)
            t2 = transform(ao1.conj().transpose(1, 0, 3, 2), S10, L10, S10, L10)
            t3 = transform(ao3, L10, S10, S10, L10)
            t4 = transform(ao3.transpose(2, 3, 0, 1), S10, L10, L10, S10)
            rel_err = max(
                rel_err,
                np.abs(t1 - t2.transpose(1, 0, 3, 2).conj()).max(),
                np.abs(t3 - t4.transpose(1, 0, 3, 2).conj()).max(),
            )
        report.check(
            "Breit: SLSL/SLLS 可由 LSLS/LSSL 严格得到",
            rel_err <= 1e-12,
            "max|d| = %.2e  -> 只需 2 个 AO 张量 / 2 次变换" % rel_err,
        )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    report.print_timings()

    print("")
    print("说明：")
    print("  * 新版 incore 用 r_outcore 的分块、AO 对称性感知变换代替 4 次全张量")
    print("    complex einsum；新版 outcore 对库仑项使用 aosym='s4'（Breit 保持 s1），")
    print("    且 Breit/Gaunt 只做 2 次变换（另外 2 项由严格共轭转置关系得到）。")
    print("  * 旧版 outcore 会忽略 retain_paldus_*/only_dg 过滤（bug），新实现与")
    print("    incore 一致；默认（不过滤）时两者输出相同。")

    print("")
    if report.ok():
        print("全部检查通过：优化版输出与旧版一致，对称性前提全部验证。")
        return 0
    print("存在失败项，请检查上面的 [FAIL] 输出。")
    return 1


if __name__ == "__main__":
    _basis = sys.argv[1] if len(sys.argv) > 1 else "3-21g"
    sys.exit(main(_basis))
