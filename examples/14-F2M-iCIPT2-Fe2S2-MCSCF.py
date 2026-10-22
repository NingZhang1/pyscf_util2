"""用 F2M_iCIPT2 作为 CI solver 计算 [2Fe-2S] 模型的 MCSCF/CASCI 能量。

FCIDUMP 来自
    "Active-space-model-for-Iron-Sulfur-Clusters"
    (20 轨道 / 30 电子, fe2s2 全 CI 参考能量 E[S=0] = -116.6056091)

F2M_iCIPT2 是 iCIPT2 家族里针对该 20 轨道 [2Fe-2S] 模型的 BUGA 版本，
输入文件格式与 iCIPT2 相同（见 iCIPT2_Loong/examples/06-F2M_iCIPT2/）。

注意: F2M_iCIPT2 本身是一个局域近似 (每个局域态最多 mxopen=5 个开壳层
电子, S 3p 外围配体只用有限局域激发), 它给出的能量高于 full CI 参考值。
本机测试结果 (mx_spintwo=10):

    cmin=1e-2, etol=1e-6, mx_ext_level=1 : -116.384058
    cmin=1e-2, etol=1e-6, mx_ext_level=2 : -116.411033
    cmin=1e-3, etol=1e-8, mx_ext_level=1 : -116.461360
    cmin=1e-3, etol=1e-8, mx_ext_level=2 : -116.524838  (本脚本默认, ~4 min)
    full CI (DMRG d=8000)                : -116.6056091

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
#   * 1e-3 : 几分钟（8 线程 + mx_ext_level=2 约 4 min），
#            RDM 与程序打印的能量一致
#   * 1e-4 : 十分钟量级；当前版本在 etol 很小时，BUGA dump 出来的 RDM
#            可能与最后打印的能量不自洽（差 ~0.3 eV），此时 solver 会
#            warning 并采用与 RDM 自洽的能量
CMIN = "1e-3"
ETOL = 1e-8

# mx_ext_level: S 3p 外围配体的局域激发等级
#   1: 官方例子设置；2: 更接近 full CI（mx_ext_level > 2 时已饱和）
MX_EXT_LEVEL = 2

# mx_spintwo: BUGA 局域态允许的最大自旋
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
        mx_ext_level=MX_EXT_LEVEL,
        mc_conv_tol=1e-7,
        mc_max_macro=32,
        run_mcscf=True,  # False: 只做一次 CASCI
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
