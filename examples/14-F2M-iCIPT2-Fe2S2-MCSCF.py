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

本脚本默认 cmin=1e-2, mode 3 + ligand level 2。

注意: 当前 F2M_iCIPT2 版本的程序接口是
    <input> <FCIDUMP> <mx_spintwo> <mode> <mx_ligand_ext_level>
            <mx_ligand_spintwo> <calculate_dm> <only_ferro_antiferro>
            <hmat_brute_force_threshold> <hmat_update> <prune_init_guess>
其中 calculate_dm 已废弃, RDM 输出由输入文件的 dumprdm 控制
(dumprdm=2 即态平均 1+2RDM, 本 driver 默认使用)。

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

# 工作目录: FCIDUMP / input / stdout / RDM csv / mo_coeff 都写在这里
WORKDIR = os.path.join(os.getcwd(), "f2m_icipt2_fe2s2")

# 局域态模式 (累积): 1/2/3
MODE = 3

# cmin: iCIPT2 的 selection threshold (输入文件里 selection=1 时生效)
#   * 1e-2 : 几十秒，粗指标
#   * 1e-3 : 几分钟，能量更低
CMIN = "1e-2"
ETOL = "1e-6"

# mode 3 才用到的配体激发参数
#   mx_ligand_ext_level: 0..6 (S3p 最多 6 个电子; 外围 S 为 0..2)
#   mx_ligand_spintwo  : 配体局域态允许的最大 2S (3 = 不限制)
MX_LIGAND_EXT_LEVEL = 2
MX_LIGAND_SPINTWO = 3

# mx_spintwo: BUGA 局域态允许的最大自旋 (None -> 模型默认值)
MX_SPINTWO = 10

# 其他新选项 (与程序命令行一一对应)
SELECTION = 1  # 1: 做 cmin selection; 0: 初猜即最终结果
ONLY_FERRO_ANTIFERRO = False
HMAT_BRUTE_FORCE_THRESHOLD = 0
HMAT_UPDATE = True
PRUNE_INIT_GUESS = True

# 每个 macro 之后对 mo_coeff 做块对角规范旋转, 使其尽量对角。
# 默认 gauge_partition="occ": 只旋转双占的 S 轨道组 (严格能量不变);
# 可选 "occ_act" 同时旋转活性 Fe 组 (只是近似不变)。
# mode 3 (LMCT) 时全部轨道都是活性轨道, 自动关闭。
GAUGE_AFTER_MACRO = True
GAUGE_PARTITION = "occ"


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
        selection=SELECTION,
        only_ferro_antiferro=ONLY_FERRO_ANTIFERRO,
        hmat_brute_force_threshold=HMAT_BRUTE_FORCE_THRESHOLD,
        hmat_update=HMAT_UPDATE,
        prune_init_guess=PRUNE_INIT_GUESS,
        gauge_after_macro=GAUGE_AFTER_MACRO,
        gauge_partition=GAUGE_PARTITION,
        mc_conv_tol=1e-7,
        mc_max_macro=8,  # 示例只跑 8 个 macro；收敛通常需要更多
        run_mcscf=True,  # False: 只做一次 CASCI
        print_summary=True,
        verbose=4,
    )

    print("")
    print("=" * 70)
    print("mode                     = %d" % MODE)
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
