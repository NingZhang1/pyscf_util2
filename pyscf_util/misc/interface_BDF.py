import pyscf
from pyscf import scf
import numpy as np
from mokit.lib.py2fch_direct import fchk
from mokit.lib.gaussian import mo_fch2py
import os
from pyscf_util.misc._parse_bdf_chkfil import (
    read_ao2somat_from_chkfil,
    ao2somat_split_based_on_irrep,
)
from pyscf_util.misc._parse_bdf_orbfile import BDFOrbParser
from pyscf_util.misc.dump_to_bdforb import dump_to_scforb


def convert_bdf_to_pyscf(
    Mol: pyscf.gto.Mole,
    mf: scf.RHF,
    chkfil_path="02S.chkfil",  # 必须保证 chkfil 文件 和 Mol 和 scforb 文件分子构型时一致的，注意，bdf 做带对称性的计算时候会调整分子构型
    scforb_path="bdf_test/02S.scforb",
    output_fch="test.fch",
    output_fch_new="test_new.fch",
    output_scforb="02S_nosymm.scforb",
    old_bdf_convention=False,
    is_casorb=False,
):
    """将 BDF 程序的分子轨道转换为 PySCF 格式的分子轨道。

    此函数实现了从 BDF 程序到 PySCF 的分子轨道转换。主要步骤包括：
    1. 读取 BDF 的对称性轨道信息
    2. 读取 BDF 的分子轨道系数
    3. 将轨道转换为 PySCF 格式
    4. 通过 MOKIT 工具进行格式转换
    5. 清理临时文件

    注意事项：
    - 输入的 chkfil 文件必须与 Mol 对象和 scforb 文件具有相同的分子构型
    - BDF 在进行对称性计算时可能会调整分子构型，需要特别注意
    - 函数会创建并删除临时的 fch 文件

    Parameters
    ----------
    Mol : pyscf.gto.Mole
        PySCF 的分子对象，包含分子构型、基组等信息
    mf : scf.RHF
        PySCF 的 RHF 计算对象
    chkfil_path : str, optional
        BDF 的 checkpoint 文件路径，默认 "02S.chkfil"
    scforb_path : str, optional
        BDF 的轨道文件路径，默认 "bdf_test/02S.scforb"
    output_fch : str, optional
        中间 fch 文件路径，默认 "test.fch"
    output_fch_new : str, optional
        最终 fch 文件路径，默认 "test_new.fch"
    output_scforb : str, optional
        输出的 BDF 轨道文件路径，默认 "02S_nosymm.scforb"
    old_bdf_convention : bool, optional
        是否使用旧的 BDF 轨道约定 (3 4 互换, 7 8 互换 (1-based))，默认 False
    is_casorb : bool, optional
        是否为 CAS 轨道，默认 False

    Returns
    -------
    numpy.ndarray
        转换后的分子轨道系数矩阵，形状为 (nao, nao)，其中 nao 为原子轨道数量

    Raises
    ------
    ValueError
        当环境变量 BDF2FCH 未设置时抛出
    """
    # 检查必要的环境变量
    BDF2FCH = os.getenv("BDF2FCH")
    if BDF2FCH is None:
        raise ValueError("BDF2FCH is not set")

    # 保存原始的 SCF 最大循环次数
    max_cycle_bak = mf.max_cycle
    # 设置 SCF 参数，但不进行实际计算
    mf.max_cycle = 1
    # 初始化分子轨道系数和能量
    mf.mo_coeff = np.zeros((Mol.nao, Mol.nao))
    mf.mo_energy = np.zeros(Mol.nao)

    # 读取 BDF 的对称性轨道信息
    ao2somat_bdf = read_ao2somat_from_chkfil(chkfil_path)
    ao2somat_bdf = ao2somat_split_based_on_irrep(ao2somat_bdf, Mol)

    # 读取并解析 BDF 轨道文件
    parser = BDFOrbParser(scforb_path)
    parser.parse_file()
    # 如果需要，转换为新的 BDF 约定
    if old_bdf_convention:
        parser.BDFold_2_new()

    # 收集轨道数据
    mo_coeffs = []
    energies = []
    occupancies = []

    # 按不可约表示处理轨道
    for irrep in range(len(Mol.irrep_name)):
        mo_coeff_tmp = ao2somat_bdf[irrep] @ parser.get_sym_data(irrep)
        mo_coeffs.append(mo_coeff_tmp)
        energies.append(parser.get_sym_energies(irrep))
        occupancies.append(parser.get_sym_occupations(irrep))

    # 合并所有不可约表示的轨道数据
    mo_coeffs = np.hstack(mo_coeffs)
    energies = np.hstack(energies)
    occupancies = np.hstack(occupancies)

    # 将轨道数据写入 BDF 格式
    dump_to_scforb(
        Mol, mo_coeffs, energies, occupancies, output_scforb, is_casorb=is_casorb
    )

    # 使用 MOKIT 工具进行格式转换
    fchk(mf, output_fch)
    os.system(f"{BDF2FCH} {output_scforb} {output_fch} {output_fch_new}")

    # 读取转换后的轨道
    mo_coeffs_bdf = mo_fch2py(output_fch_new)

    # 恢复原始的 SCF 最大循环次数
    mf.max_cycle = max_cycle_bak

    # 清理临时文件
    os.system(f"rm {output_fch}")
    os.system(f"rm {output_fch_new}")

    return mo_coeffs_bdf


if __name__ == "__main__":
    # 示例：使用 C10H8 分子进行测试
    GEOMETRY = """
C              4.598900802567      -0.000000000000      -1.338594114377
C              4.598900802567      -0.000000000000       1.338594114377
H              6.383749016945      -0.000000000000       2.354252931327
H              6.383749016945      -0.000000000000      -2.354252931327
C              2.352149666151      -0.000000000000      -2.650354117785
C              2.352149666151      -0.000000000000       2.650354117785
C             -0.000000000000      -0.000000000000      -1.355251105302
C              0.000000000000       0.000000000000       1.355251105302
H              2.348156674849      -0.000000000000      -4.705947173482
H              2.348156674849      -0.000000000000       4.705947173482
C             -2.352149666151       0.000000000000      -2.650354117785
C             -2.352149666151       0.000000000000       2.650354117785
C             -4.598900802567       0.000000000000      -1.338594114377
C             -4.598900802567       0.000000000000       1.338594114377
H             -2.348156674849       0.000000000000      -4.705947173482
H             -2.348156674849       0.000000000000       4.705947173482
H             -6.383749016945       0.000000000000      -2.354252931327
H             -6.383749016945       0.000000000000       2.354252931327
"""
    # 构建分子对象
    Mol = pyscf.gto.Mole()
    Mol.atom = GEOMETRY
    Mol.basis = "cc-pvdz"
    Mol.symmetry = "d2h"
    Mol.spin = 0
    Mol.charge = 0
    Mol.verbose = 4
    Mol.unit = "bohr"
    Mol.build()

    # 构建 RHF 计算对象
    mf = scf.RHF(Mol)

    # 转换为 BDF 格式并获取分子轨道
    mo_coeffs_bdf = convert_bdf_to_pyscf(Mol, mf)

    # 检查轨道的正交性
    ovlp = Mol.intor("int1e_ovlp")
    ovlp_mo = mo_coeffs_bdf.T @ ovlp @ mo_coeffs_bdf

    # 打印重叠矩阵的对角元素（应该接近 1.0）
    print(np.diag(ovlp_mo))

    # 可选：使用转换后的轨道进行 SCF 计算
    # dm_init = mf.make_rdm1(mo_coeffs_bdf)
    # mf.kernel(dm0=dm_init)  # 应该在一两个循环内收敛
