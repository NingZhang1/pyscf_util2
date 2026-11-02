"""用 F2M_iCIPT2 作为 CI solver 计算 [2Fe-2S] 模型的 MCSCF/CASCI 能量。

FCIDUMP 来自
    "Active-space-model-for-Iron-Sulfur-Clusters"
    (20 轨道 / 30 电子, fe2s2 全 CI 参考能量 E[S=0] = -116.6056091)

F2M_iCIPT2 是 iCIPT2 家族里针对该 20 轨道 [2Fe-2S] 模型的 BUGA 版本，
输入文件格式与 iCIPT2 相同（见 iCIPT2_Loong/examples/06-F2M_iCIPT2/）。

局域态模式 (与 Fe4S4 版本一致, 累积):

    mode 1 : 只有自旋激发            Fe(2+) / Fe(3+)
    mode 2 : + MMCT                  Fe(+) / Fe(4+) (d^4 high spin)
    mode 3 : + LMCT                  配体 S 3p 局域激发

注意: F2M_iCIPT2 是一个局域近似 (每个局域态只保留最高自旋多重态,
mxopen=5), 能量高于 full CI 参考值。本机测试结果 (mx_spintwo=10,
cmin=1e-2, etol=1e-6, 纯 CASCI):

    mode 1                          : -115.953720  (Fe2S2 中两个 Fe 都是
                                      Fe(3+), 只有一个高自旋组态)
    mode 2                          : -115.955770
    mode 3, ligand level 1          : -116.383987
    mode 3, ligand level 2 或 3     : -116.410994
    full CI (DMRG d=8000)           : -116.6056091

cmin=1e-3 时 mode 3 + ligand level 2 给出 -116.524975 (ncsf=6388)。
本脚本默认 cmin=1e-2, mode 3 + ligand level 2。

运行前设置可执行文件路径:
    export F2M_ICIPT2=/path/to/iCIPT2_Loong/bin/F2M_iCIPT2_Fe2S2.exe

用法:
    python 14-F2M-iCIPT2-Fe2S2-MCSCF.py
"""

import os
import time

import numpy

from pyscf_util.MeanField import f2micipt2

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
        "Fe2S2",
        "fe2s2",
    )
)
FCIDUMP = os.environ.get("FE2S2_FCIDUMP", _DEFAULT_FCIDUMP)

# ------------------------------------------------------------------
# 计算设置
# ------------------------------------------------------------------

# F2M_iCIPT2 可执行文件 (由环境变量 F2M_ICIPT2 或 PATH 提供)
EXECUTABLE = os.environ.get("F2M_ICIPT2")

# 工作目录: FCIDUMP / input / stdout / RDM csv 都写在这里
WORKDIR = os.path.join(os.getcwd(), "f2m_icipt2_fe2s2")

# cmin: iCIPT2 的 selection threshold。
#   * 1e-2 : 几十秒，粗指标
#   * 1e-3 : 几分钟，RDM 与程序打印的能量一致
#   * 1e-4 : 十分钟量级；当前版本在 etol 很小时，BUGA dump 出来的 RDM
#            可能与最后打印的能量不自洽（差 ~0.3 eV），此时 solver 会
#            warning 并采用与 RDM 自洽的能量
CMIN = "1e-2"
ETOL = "1e-6"

# 局域态模式 (与 Fe4S4 一致):
#   1: 只有自旋激发 (Fe2+/Fe3+)
#   2: + MMCT (Fe+/Fe4+)
#   3: + LMCT (配体激发)
MODE = 3

# mode 3 才用到的配体激发参数
#   mx_ligand_ext_level: 0..6 (S3p 最多 6 个电子; 外围 S 为 0..2)
#   mx_ligand_spintwo  : 配体局域态允许的最大 2S (3 才允许 p^3 的 quartet)
MX_LIGAND_EXT_LEVEL = 2
MX_LIGAND_SPINTWO = 2

# mx_spintwo: BUGA 局域态允许的最大自旋 (None -> 模型默认值)
MX_SPINTWO = 10


def main():
    t0 = time.time()
    mc = f2micipt2.kernel(
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
        mc_max_macro=8,  # 示例只跑 8 个 macro；收敛通常需要更多
        run_mcscf=True,  # False: 只做一次 CASCI
        print_summary=True,
        verbose=4,
    )

    print("")
    print("=" * 70)
    print("F2M_iCIPT2 / FCI energy  = %.10f" % mc.e_tot)
    print("ncas / nelecas           = %d / %s" % (mc.ncas, mc.nelecas))
    print("reference (DMRG d=8000)  = -116.6056091")
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
