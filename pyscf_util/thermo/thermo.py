import numpy as np
from pyscf import gto, scf, dft, hessian
from pyscf.hessian import thermo

import numpy
from pyscf.data import nist


def compute_moments_of_inertia(mol):
    """
    计算分子主转动惯量，单位 amu * Bohr^2。
    返回按升序排列的三个主转动惯量。
    != rotation constant !
    """
    coords = mol.atom_coords()  # (natm, 3) 单位 Bohr
    masses = np.array(mol.atom_mass_list(isotope_avg=True))  # 单位 amu

    # 1. 计算质心
    total_mass = masses.sum()
    com = (masses[:, None] * coords).sum(axis=0) / total_mass

    # 2. 相对质心坐标
    r = coords - com

    # 3. 构建惯性张量
    I = np.zeros((3, 3))
    for i in range(mol.natm):
        m = masses[i]
        x, y, z = r[i]
        I[0, 0] += m * (y * y + z * z)
        I[1, 1] += m * (x * x + z * z)
        I[2, 2] += m * (x * x + y * y)
        I[0, 1] -= m * x * y
        I[0, 2] -= m * x * z
        I[1, 2] -= m * y * z
    I[1, 0] = I[0, 1]
    I[2, 0] = I[0, 2]
    I[2, 1] = I[1, 2]

    # 4. 对角化，得到主转动惯量
    eigvals = np.linalg.eigvalsh(I)
    return np.sort(eigvals)


def generate_kinp(
    mol,
    mf,
    hessian_matrix=None,
    symmetry_number=1,
    electronic_degeneracy=1,
    output_file="output.kinp",
):
    """
    根据 PySCF 计算结果生成 KiSThelP 的 .kinp 输入文件。

    参数:
        mol: pyscf.gto.Mole 对象
        mf: 平均场对象 (SCF/DFT)，需已收敛，包含 e_tot
        hessian_matrix: 可选，Hessian 矩阵。若为 None，则自动计算。
        symmetry_number: 转动对称数 (整数)
        electronic_degeneracy: 电子态简并度 (整数)
        output_file: 输出文件名

    返回:
        写入的文件路径
    """
    # 1. 计算 Hessian（如果需要）
    if hessian_matrix is None:
        # 根据 mf 类型选择 Hessian 计算方式
        if isinstance(mf, dft.RKS) or isinstance(mf, dft.UKS):
            hessian_matrix = hessian.RKS(mf).kernel()
        elif isinstance(mf, scf.RHF) or isinstance(mf, scf.UHF):
            hessian_matrix = hessian.RHF(mf).kernel()
        else:
            raise TypeError("不支持的 mf 类型，请手动提供 hessian_matrix")

    # 2. 频率分析
    freq_results = thermo.harmonic_analysis(mol, hessian_matrix)
    freqs = freq_results["freq_wavenumber"]  # 包含平动转动

    natm = mol.natm
    # 判断线性：线性分子有 3N-5 个振动模式，非线性 3N-6
    # 总频率数应为 3N（包括平动转动），振动模式数 = len(freqs) - 6 或 -5
    # 实际上 harmonic_analysis 返回的 freqs 长度是 3N，前6或5个接近0
    # 我们可以通过检查最小的几个频率是否接近0来判断，但更简单：根据分子几何判断线性
    # 这里通过计算转动惯量：线性分子有一个主转动惯量为0
    moments = compute_moments_of_inertia(mol)
    # 如果最小转动惯量小于某个阈值，认为是线性
    if moments[0] < 1e-6:
        is_linear = True
        n_vib = 3 * natm - 5
    else:
        is_linear = False
        n_vib = 3 * natm - 6

    # 提取振动频率：通常 freqs 已经排序，前几个是平动转动，接近0
    # 但为了保险，我们取绝对值最大的 n_vib 个？或者直接取最后 n_vib 个
    # harmonic_analysis 返回的 freqs 是排序的，但可能包含负值（虚频）
    # 我们取绝对值大于某个阈值的？简单：取最后 n_vib 个
    vib_freqs = freqs[-n_vib:] if n_vib > 0 else []

    # 3. 总质量
    total_mass = sum(mol.atom_mass_list(isotope_avg=True))

    # 4. 电子能量
    potential_energy = mf.e_tot

    # 5. 转动惯量（非零的）
    if is_linear:
        # 线性分子，KiSThelP 可能只需要一个转动惯量？查看文档：对于线性分子，转动惯量只有一个非零值
        # 通常线性分子有两个相等的转动惯量，另一个为0。KiSThelP 可能要求提供非零的那个？或者全部？文档说 *MOMENT OF INERTIA 列表
        # 对于线性分子，可能只需要一个值。但为安全，提供非零值。
        nonzero_moments = moments[moments > 1e-6]
        moments_of_inertia = nonzero_moments
    else:
        moments_of_inertia = moments  # 三个主转动惯量

    # 6. 构建 .kinp 内容
    lines = []
    lines.append("*MASS")
    lines.append(f"{total_mass:.6f}")
    lines.append("*END")
    lines.append("")
    lines.append("*FREQUENCIES")
    for f in vib_freqs:
        lines.append(f"{f:.4f}")
    lines.append("*END")
    lines.append("")
    lines.append("*LINEAR")
    lines.append("linear" if is_linear else "not linear")
    lines.append("*END")
    lines.append("")
    lines.append("*NUMBER OF SYMMETRY")
    lines.append(str(symmetry_number))
    lines.append("*END")
    lines.append("")
    lines.append("*MOMENT OF INERTIA")
    for mi in moments_of_inertia:
        lines.append(f"{mi:.6f}")
    lines.append("*END")
    lines.append("")
    lines.append("*POTENTIAL ENERGY")
    lines.append(f"{potential_energy:.10f}")
    lines.append("*END")
    lines.append("")
    lines.append("*ELECTRONIC DEGENERACY")
    lines.append(str(electronic_degeneracy))
    lines.append("*END")

    content = "\n".join(lines) + "\n"

    with open(output_file, "w") as f:
        f.write(content)

    print(f"✅ 已生成 {output_file}")
    return output_file


if __name__ == "__main__":

    from numerical_hessian import *
    from pyscf.geomopt.geometric_solver import optimize

    # 1. 定义分子
    mol = gto.M(
        atom="""
        O  0.000000  0.000000  0.000000
        H  0.000000  0.000000  0.9572
        H  0.926627  0.000000 -0.239987
        """,
        basis="6-31g(d,p)",
        unit="Angstrom",
        verbose=0,
    )

    mf = scf.RHF(mol)

    mol = optimize(
        mf,
        maxsteps=100,
        conv_params={
            "convergence_energy": 1e-6,
            "convergence_grms": 3e-4,
            "convergence_gmax": 4.5e-4,
            "convergence_drms": 1.2e-3,
            "convergence_dmax": 1.8e-3,
        },
    )

    # 初始核坐标，展平为一维向量，单位 Bohr
    x0 = mol.atom_coords(unit="Bohr").ravel()
    natm = mol.natm

    # 原子质量 (amu)。也可以用 mol.atom_mass_list(isotope_avg=True)
    masses = np.array(mol.atom_mass_list(isotope_avg=True))

    # 2. 定义能量函数：输入核坐标向量 x，返回 RHF 能量
    def energy(x):
        # 更新分子几何
        mol.set_geom_(x.reshape(-1, 3), unit="Bohr")
        # 每次都重新做 RHF，保证从当前几何出发
        mf = scf.RHF(mol)
        mf.conv_tol = 1e-14
        mf.conv_tol_grad = 1e-10
        mf.max_cycle = 200
        e = mf.kernel()
        return e

    # ---------- 解析 hess ----------

    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-10
    mf.max_cycle = 200
    mf.kernel()

    # 解析 Hessian：单位 Hartree/Bohr^2，shape (natm, natm, 3, 3)

    hess_ana = mf.Hessian().kernel()

    # ---------- numerical hess ----------

    grad, hess = numerical_grad_hess(energy, x0, step=1e-3, verbose=True)
    hess = hess.reshape(natm, 3, natm, 3).transpose(0, 2, 1, 3)

    print(grad)

    freq_results_num = thermo.harmonic_analysis(mol, hess)
    freq_results_ana = thermo.harmonic_analysis(mol, hess_ana)

    # print(freq_results_ana)
    # print(freq_results_num)

    moments = compute_moments_of_inertia(mol)

    print(moments)
    # print(thermo.rotation_const(mol.atom_mass_list(isotope_avg=True), mol.atom_coords(), unit="GHz"))

    generate_kinp(
        mol,
        mf,
        hessian_matrix=hess,
        symmetry_number=thermo.rotational_symmetry_number(mol),
        electronic_degeneracy=1,
        output_file="output.kinp",
    )
