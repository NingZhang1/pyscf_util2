#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""四分量积分：ASCII head + binary 积分（供 Loong iCIPT2 使用）。

用法::

    python examples/18-4C-ints-binary.py                # 默认 F/sto-3g（秒级）
    python examples/18-4C-ints-binary.py unc-cc-pvdz-dk # 更大基组（分钟级）

产物（以 ``FCIDUMP_F`` 为例）：

* ``FCIDUMP_F``      : 只含 ASCII head（NORB/NELEC/MS2/ORBSYM/ISYM）
* ``FCIDUMP_F.bin``  : 全部积分（2e + 1e + ECORE），每条 16 字节（实数）

C++ 端（``iCIPT2_Breit_D2H_ERICompressed``）用第 2 个参数读 head、第 3 个
参数读 binary 积分::

    ./iCIPT2_Breit_D2H_ERI_compressed.exe input FCIDUMP_F FCIDUMP_F.bin

脚本会：

1. 用 ``FCIDUMP_Rela4C_SU2_opt`` 生成传统 ASCII FCIDUMP（参考）；
2. 用 ``FCIDUMP_Rela4C_SU2_binary`` 生成 head + binary；
3. 逐条比较 ASCII 与 binary（记录集合、积分值、head）；
4. 演示 ``binary_to_fcidump`` / ``fcidump_to_binary`` 双向转换。

任何检查失败时以退出码 1 结束。
"""

import copy
import os
import shutil
import sys
import tempfile

import numpy as np

from pyscf import gto, scf

from pyscf_util.Relativisitc.integral_4C_binary import (
    FCIDUMP_Rela4C_SU2_binary,
    binary_summary,
    binary_to_fcidump,
    compare_binary_with_fcidump,
    fcidump_to_binary,
)
from pyscf_util.Relativisitc.integral_4C_opt import FCIDUMP_Rela4C_SU2_opt

ATOM = "F 0 0 0"  # D2h 原子，对应 *_D2H 程序


def main(basis="3-21g"):
    workdir = tempfile.mkdtemp(prefix="rela4c_binary_example_")
    print("workdir = %s" % workdir)
    try:
        mol = gto.M(
            atom=ATOM,
            basis=basis,
            verbose=0,
            charge=-1,
            spin=0,
            symmetry="d2h",
        )
        mol.build()
        mf = scf.dhf.RDHF(mol)
        mf.conv_tol = 1e-10
        mf.kernel()
        mf.with_breit = True
        mo0 = np.array(mf.mo_coeff)

        # ---------------- 1. ASCII reference (what the codes read today) ---
        ref = os.path.join(workdir, "FCIDUMP_ascii")
        mf_ref = copy.copy(mf)
        mf_ref.mo_coeff = mo0.copy()
        FCIDUMP_Rela4C_SU2_opt(mol, mf_ref, True, filename=ref, mode="incore")
        print(
            "[1] ASCII reference : %s (%.1f MB)" % (ref, os.path.getsize(ref) / 2.0**20)
        )

        # ---------------- 2. head + binary --------------------------------
        head = os.path.join(workdir, "FCIDUMP_F")
        mf_bin = copy.copy(mf)
        mf_bin.mo_coeff = mo0.copy()
        FCIDUMP_Rela4C_SU2_binary(
            mol,
            mf_bin,
            True,
            filename=head,
            mode="incore",
            # write_ascii_integrals=True,  # also keep the classic ASCII records
        )
        binary = head + ".bin"
        summary = binary_summary(binary, is_complex=False)
        print(
            "[2] binary          : %s (%.1f MB)"
            % (binary, os.path.getsize(binary) / 2.0**20)
        )
        print("    records         : %s" % summary)
        print(
            "    size ratio      : %.2f%%"
            % (100.0 * os.path.getsize(binary) / os.path.getsize(ref))
        )

        # ---------------- 3. binary == ASCII ------------------------------
        result = compare_binary_with_fcidump(
            binary, ref, is_complex=False, tol=1e-10, verbose=True
        )
        print("[3] binary == ASCII : OK (n_2e = %d)" % result["n_2e"])

        # ---------------- 4. conversions ----------------------------------
        back = os.path.join(workdir, "back_to_ascii")
        binary_to_fcidump(binary, back, head_filename=head)
        split = os.path.join(workdir, "resplit.bin")
        stats = fcidump_to_binary(ref, head + ".resplit", split)
        print("[4] binary -> ASCII : %s" % back)
        print("    ASCII  -> binary : %s (%s)" % (split, stats))

        print(
            "\nNow run, e.g.:\n"
            "  ./iCIPT2_Breit_D2H_ERI_compressed.exe input %s %s\n"
            "  ./iCIPT2_Breit_D2H_ERI_compressed.exe input %s        # ASCII reference"
            % (os.path.basename(head), os.path.basename(binary), os.path.basename(ref))
        )
        return result
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    try:
        main(sys.argv[1] if len(sys.argv) > 1 else "3-21g")
    except AssertionError as exc:
        print("FAILED:", exc)
        sys.exit(1)
