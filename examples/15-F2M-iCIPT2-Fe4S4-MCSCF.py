"""用 F2M_iCIPT2_Fe4S4 作为 CI solver 计算 [4Fe-4S] 模型的 MCSCF/CASCI 能量。

FCIDUMP 来自
    "Active-space-model-for-Iron-Sulfur-Clusters"
    (36 轨道 / 54 电子, Fe4S4 全 CI 参考能量 E[S=0] = -327.2396369)

F2M_iCIPT2_Fe4S4 是 F2M_iCIPT2_Fe2S2 的 [4Fe-4S] 版本，支持三种局域态
模式 (见 pyscf_util/MeanField/f2micipt2_fe4s4.py)：

    mode 1 : 只有自旋激发            Fe(2+) / Fe(3+)
    mode 2 : + MMCT                  Fe(+) / Fe(4+) (d^4 high spin)
    mode 3 : + LMCT                  配体 S 3p 局域激发，
                                     重数由 MX_LIGAND_SPINTWO 限制

注意: F2M_iCIPT2 是一个局域近似，能量高于 full CI 参考值；cmin 越小、
模式越"大"，能量越低 (但代价越高)。

本机测试结果 (mx_spintwo=20, cmin=1e-1, etol=1e-4; 只用于验证程序能跑通):

    mode 1 (spin only)   : -326.215820   (~2 min)
    mode 2 (MMCT)        : -326.215820   (cmin=1e-1 时 MMCT 态未被选到)
    mode 3 (LMCT)        : -326.414887
    full CI (DMRG d=4000): -327.2396369

运行前设置可执行文件路径:
    export F2M_ICIPT2_FE4S4=/path/to/iCIPT2_Loong/bin/F2M_iCIPT2_Fe4S4.exe

用法:
    python 15-F2M-iCIPT2-Fe4S4-MCSCF.py
"""

import os
import time

import numpy

from pyscf_util.MeanField import f2micipt2_fe4s4

# ------------------------------------------------------------------
# FCIDUMP
# ------------------------------------------------------------------

_here = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_FCIDUMP = os.path.abspath(
    os.path.join(
        _here,
        "..",
        "..",
        "Active-space-model-for-Iron-Sulfur-Clusters",
        "Fe2S2_and_Fe4S4",
        "Fe4S4",
        "fe4s4",
    )
)
FCIDUMP = os.environ.get("FE4S4_FCIDUMP", _DEFAULT_FCIDUMP)

# ------------------------------------------------------------------
# 计算设置
# ------------------------------------------------------------------

# F2M_iCIPT2_Fe4S4 可执行文件 (由 F2M_ICIPT2_FE4S4 环境变量或 PATH 提供)
EXECUTABLE = os.environ.get("F2M_ICIPT2_FE4S4")

# 工作目录
WORKDIR = os.path.join(os.getcwd(), "f2m_icipt2_fe4s4")

# 局域态模式: 1 / 2 / 3
MODE = 1

# cmin: iCIPT2 的 selection threshold
#   * 1e-1 : 一分钟量级，只用于验证程序能跑通
#   * 1e-2 : 更贵；会得到更低的能量
CMIN = "1e-1"
ETOL = 1e-4

# mx_spintwo: BUGA 局域态允许的最大自旋
#   [4Fe-4S] 有 4 个 Fe，中间自旋最大可以到 4*5/2*2 = 20
MX_SPINTWO = 20

# mode 3 才用到的配体激发参数
MX_LIGAND_EXT_LEVEL = 1
MX_LIGAND_SPINTWO = 1


def main():
    t0 = time.time()
    mc = f2micipt2_fe4s4.kernel(
        FCIDUMP,
        executable=EXECUTABLE,
        workdir=WORKDIR,
        spin=0,  # singlet
        cmin=CMIN,
        etol=ETOL,
        mx_spintwo=MX_SPINTWO,
        mode=MODE,
        mx_ligand_ext_level=MX_LIGAND_EXT_LEVEL,
        mx_ligand_spintwo=MX_LIGAND_SPINTWO,
        mc_conv_tol=1e-7,
        mc_max_macro=32,
        run_mcscf=True,  # False: 只做一次 CASCI
        verbose=4,
    )

    print("")
    print("=" * 70)
    print("mode                     = %d" % MODE)
    print("F2M_iCIPT2_Fe4S4 energy  = %.10f" % mc.e_tot)
    print("ncas / nelecas           = %d / %s" % (mc.ncas, mc.nelecas))
    print("reference (DMRG d=4000)  = -327.2396369")
    print("wall time                = %.1f s" % (time.time() - t0))
    print("=" * 70)

    # 自然轨道占据数
    dm1 = mc.fcisolver.make_rdm1()
    occ = numpy.linalg.eigvalsh(dm1)[::-1]
    print("natural occupations:")
    print(numpy.array2string(occ, precision=6, suppress_small=True))
    return mc


if __name__ == "__main__":
    main()
