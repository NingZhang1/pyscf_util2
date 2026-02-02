import numpy as np
from pyscf import lo, gto, scf
from pyscf.lo import boys
import warnings
import pyscf


def random_orthogonal_matrix(n):
    """
    生成n×n的随机正交矩阵（实数）
    """
    # 生成随机矩阵
    A = np.random.randn(n, n)
    # QR分解得到正交矩阵
    Q, R = np.linalg.qr(A)
    return Q


def atomic_canonicalization(
    mol, fock, mo_coeff, SCF=None, threshold1=0.7, threshold2=0.99
):
    """
    执行原子正则化的函数

    参数:
    mol: PySCF分子对象
    fock: Fock矩阵 (AO基)
    ovlp: 重叠矩阵 (AO基)
    mo_coeff: 输入分子轨道系数
    threshold1: 第一个阈值，用于判断轨道局域化程度
    threshold2: 第二个阈值，用于判断轨道局域化程度

    返回:
    mo_atm_can: 原子正则化轨道系数，如果失败返回False
    """

    # (0) apply random matrix

    # 判断轨道对称性

    if mol.symmetry == False:
        norb = mo_coeff.shape[1]
        R = random_orthogonal_matrix(norb)
        mo_processed = mo_coeff @ R
    else:
        OrbSym = pyscf.symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, mo_coeff)
        unique_sym = list(set(OrbSym))
        mo_processed = mo_coeff.copy()

        for sym in unique_sym:
            indx = [np.where(np.array(OrbSym) == sym)[0]]
            n_orb = len(indx[0])
            if n_orb > 1:
                R = random_orthogonal_matrix(n_orb)
                mo_processed[:, indx[0]] = mo_processed[:, indx[0]] @ R

    # (1) Boys局域化
    print("进行Boys局域化...")
    if SCF is None:
        mo_boys = boys.Boys(mol, mo_processed).kernel()
    else:
        # import pyscf
        pm = pyscf.lo.PM(mol, mo_processed, SCF)
        mo_boys = pm.kernel()

    # (2) 计算每个局域轨道的<r>并分析
    print("分析轨道局域化程度...")
    atom_coords = mol.atom_coords()
    nocc = mo_boys.shape[1]

    # 计算轨道中心 <r> = <phi|r|phi>
    orb_centers = np.zeros((nocc, 3))

    rx, ry, rz = mol.intor_symmetric("int1e_r")
    ints = {
        "x": rx,
        "y": ry,
        "z": rz,
    }

    for i in range(nocc):
        # 轨道密度矩阵: |mo_i><mo_i|
        dm_i = np.outer(mo_boys[:, i], mo_boys[:, i])
        # 计算 <x>, <y>, <z>
        for j, coord in enumerate(["x", "y", "z"]):
            # ints = mol.intor_symmetric(f"int1e_r{coord}", hermi=1)
            orb_centers[i, j] = np.trace(dm_i @ ints[coord])

    # 找出每个轨道最近的两个原子及其距离
    nearest_atoms = []
    distance_ratios = []
    orbital_warnings = []

    for i in range(nocc):
        distances = []
        for atom_idx in range(len(atom_coords)):
            dist = np.linalg.norm(orb_centers[i] - atom_coords[atom_idx])
            distances.append((dist, atom_idx))

        # 按距离排序
        distances.sort(key=lambda x: x[0])
        # print(distances)
        nearest_dist, nearest_atom = distances[0]
        second_dist, second_atom = distances[1]

        nearest_atoms.append(nearest_atom)

        ratio = nearest_dist / second_dist

        # if nearest_dist > 0:
        #     ratio = second_dist / nearest_dist
        # else:
        #     ratio = float("inf")

        distance_ratios.append(ratio)

        # 检查局域化程度
        if threshold1 < ratio <= threshold2:
            warning_msg = f"轨道 {i} 难以区分: 最近原子 {nearest_atom} 与次近原子 {second_atom} 的距离比 = {ratio:.3f}"
            warnings.warn(warning_msg)
            orbital_warnings.append(warning_msg)
        elif threshold2 < ratio <= 1.0:
            warning_msg = f"轨道 {i} 无法局域化: 最近原子 {nearest_atom} 与次近原子 {second_atom} 的距离比 = {ratio:.3f}"
            warnings.warn(warning_msg)
            orbital_warnings.append(warning_msg)
            return False, None, orbital_warnings
        else:
            if ratio > 1.0:
                raise RuntimeError

    # (3) 按所属原子分组
    print("按原子分组轨道...")
    atom_groups = {}
    for i, atom_idx in enumerate(nearest_atoms):
        if atom_idx not in atom_groups:
            atom_groups[atom_idx] = []
        atom_groups[atom_idx].append(i)

    # sort the atom groups by atom index

    atom_groups = dict(sorted(atom_groups.items()))

    norb_atm = []

    # 打印分组信息

    for atom_idx, orb_indices in atom_groups.items():
        print(f"原子 {atom_idx}: 轨道 {orb_indices}")
        norb_atm.append(len(orb_indices))

    # (4) 计算每个原子组的Fock投影并对角化
    print("计算原子正则化轨道...")
    mo_atm_can = np.zeros_like(mo_boys)

    # 将Fock矩阵变换到Boys轨道基
    fock_boys = mo_boys.T @ fock @ mo_boys

    current_col = 0
    for atom_idx, orb_indices in atom_groups.items():
        n_orbs = len(orb_indices)
        if n_orbs == 0:
            continue

        # 提取该原子组的Fock子矩阵
        indices = np.array(orb_indices)
        fock_sub = fock_boys[np.ix_(indices, indices)]

        # 对角化Fock子矩阵
        eigvals, eigvecs = np.linalg.eigh(fock_sub)

        print("原子 {} 的轨道能级: {}".format(atom_idx, eigvals))

        # 将原子正则化轨道放回正确位置
        mo_atm_can[:, current_col : current_col + n_orbs] = (
            mo_boys[:, indices] @ eigvecs
        )
        current_col += n_orbs

    # (5) 返回原子正则化轨道
    print("原子正则化完成!")
    return mo_atm_can, norb_atm, orbital_warnings


def dimer_canonicalization(
    mol, fock, mo_coeff, SCF=None, threshold1=0.7, threshold2=0.99
):
    """
    执行原子正则化的函数

    参数:
    mol: PySCF分子对象
    fock: Fock矩阵 (AO基)
    ovlp: 重叠矩阵 (AO基)
    mo_coeff: 输入分子轨道系数
    threshold1: 第一个阈值，用于判断轨道局域化程度
    threshold2: 第二个阈值，用于判断轨道局域化程度

    返回:
    mo_atm_can: 原子正则化轨道系数，如果失败返回False
    """

    # (0) apply random matrix

    # 判断轨道对称性

    if mol.symmetry == False:
        norb = mo_coeff.shape[1]
        R = random_orthogonal_matrix(norb)
        mo_processed = mo_coeff @ R
    else:
        OrbSym = pyscf.symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, mo_coeff)
        unique_sym = list(set(OrbSym))
        mo_processed = mo_coeff.copy()

        for sym in unique_sym:
            indx = [np.where(np.array(OrbSym) == sym)[0]]
            n_orb = len(indx[0])
            if n_orb > 1:
                R = random_orthogonal_matrix(n_orb)
                mo_processed[:, indx[0]] = mo_processed[:, indx[0]] @ R

    # (1) Boys局域化

    print("进行Boys局域化...")

    if SCF is None:
        mo_boys = boys.Boys(mol, mo_processed).kernel()
    else:
        # import pyscf
        pm = pyscf.lo.PM(mol, mo_processed, SCF)
        pm.max_iters = 1024
        mo_boys = pm.kernel()

    # (2) 计算每个局域轨道的<r>并分析

    print("分析轨道局域化程度...")

    atom_coords = mol.atom_coords()
    nocc = mo_boys.shape[1]

    # 计算轨道中心 <r> = <phi|r|phi>

    orb_centers = np.zeros((nocc, 3))

    rx, ry, rz = mol.intor_symmetric("int1e_r")
    ints = {
        "x": rx,
        "y": ry,
        "z": rz,
    }

    for i in range(nocc):
        # 轨道密度矩阵: |mo_i><mo_i|
        dm_i = np.outer(mo_boys[:, i], mo_boys[:, i])
        # 计算 <x>, <y>, <z>
        for j, coord in enumerate(["x", "y", "z"]):
            # ints = mol.intor_symmetric(f"int1e_r{coord}", hermi=1)
            orb_centers[i, j] = np.trace(dm_i @ ints[coord])

    # 找出每个轨道最近的两个原子及其距离

    nearest_atoms = []
    distance_ratios = []
    orbital_warnings = []

    for i in range(nocc):
        distances = []
        for atom_idx in range(len(atom_coords)):
            dist = np.linalg.norm(orb_centers[i] - atom_coords[atom_idx])
            distances.append((dist, atom_idx))

        # 按距离排序
        distances.sort(key=lambda x: x[0])
        # print(distances)
        nearest_dist, nearest_atom = distances[0]
        second_dist, second_atom = distances[1]
        third_dist, third_atom = distances[2]

        nearest_atoms.append(nearest_atom)

        if (nearest_atom // 2) != (second_atom // 2):
            print(
                f"轨道 {i} 最近的两个原子 {nearest_atom} 和 {second_atom} 不属于同一个dimer!"
            )
            return False, None, orbital_warnings

        # ratio = nearest_dist / second_dist
        ratio = second_dist / third_dist

        # if nearest_dist > 0:
        #     ratio = second_dist / nearest_dist
        # else:
        #     ratio = float("inf")

        distance_ratios.append(ratio)

        # 检查局域化程度
        if threshold1 < ratio <= threshold2:
            warning_msg = f"轨道 {i} 难以区分: 次原子 {second_atom} 与次次近原子 {third_atom} 的距离比 = {ratio:.3f}"
            warnings.warn(warning_msg)
            orbital_warnings.append(warning_msg)
        elif threshold2 < ratio <= 1.0:
            warning_msg = f"轨道 {i} 难以区分: 次近原子 {second_atom} 与次次近原子 {third_atom} 的距离比 = {ratio:.3f}"
            warnings.warn(warning_msg)
            orbital_warnings.append(warning_msg)
            return False, None, orbital_warnings
        else:
            if ratio > 1.0:
                raise RuntimeError

    # (3) 按所属原子分组
    print("按dimer分组轨道...")
    atom_groups = {}
    for i, atom_idx in enumerate(nearest_atoms):
        dimer_idx = atom_idx // 2
        if dimer_idx not in atom_groups:
            atom_groups[dimer_idx] = []
        atom_groups[dimer_idx].append(i)

    # sort the atom groups by atom index

    atom_groups = dict(sorted(atom_groups.items()))

    norb_atm = []

    # 打印分组信息

    for atom_idx, orb_indices in atom_groups.items():
        print(f"dimer {atom_idx}: 轨道 {orb_indices}")
        norb_atm.append(len(orb_indices))

    # (4) 计算每个原子组的Fock投影并对角化
    print("计算原子正则化轨道...")
    mo_atm_can = np.zeros_like(mo_boys)

    # 将Fock矩阵变换到Boys轨道基
    fock_boys = mo_boys.T @ fock @ mo_boys

    current_col = 0
    for atom_idx, orb_indices in atom_groups.items():
        n_orbs = len(orb_indices)
        if n_orbs == 0:
            continue

        # 提取该原子组的Fock子矩阵
        indices = np.array(orb_indices)
        fock_sub = fock_boys[np.ix_(indices, indices)]

        # 对角化Fock子矩阵
        eigvals, eigvecs = np.linalg.eigh(fock_sub)

        print("原子 {} 的轨道能级: {}".format(atom_idx, eigvals))

        # 将原子正则化轨道放回正确位置
        mo_atm_can[:, current_col : current_col + n_orbs] = (
            mo_boys[:, indices] @ eigvecs
        )
        current_col += n_orbs

    # (5) 返回原子正则化轨道
    print("原子正则化完成!")
    return mo_atm_can, norb_atm, orbital_warnings


# 使用示例
def example_usage():
    # 创建示例分子
    mol = gto.Mole()
    mol.atom = """
    He 0.0 0.0 0.0
    He 0.0 0.0 3.0
    He 0.0 0.0 6.0
    He 0.0 0.0 9.0
    He 0.0 0.0 12.0
    He 0.0 0.0 15.0
    """
    mol.basis = "cc-pvdz"
    mol.symmetry = False
    mol.verbose = 4
    mol.build()

    # 进行HF计算获取Fock矩阵和轨道
    mf = scf.RHF(mol)
    mf.kernel()

    # 获取Fock矩阵和重叠矩阵
    fock = mf.get_fock()
    # ovlp = mf.get_ovlp()
    mo_coeff = mf.mo_coeff

    # 执行原子正则化

    mo_atm_can, warnings = atomic_canonicalization(mol, fock, mo_coeff[:, :6])
    mo_atm_can, warnings = atomic_canonicalization(mol, fock, mo_coeff[:, 6:])

    if mo_atm_can is not False:
        print("原子正则化轨道形状:", mo_atm_can.shape)
        # 可以继续使用原子正则化轨道进行后续计算
    else:
        print("原子正则化失败")
        print("警告信息:", warnings)


if __name__ == "__main__":
    example_usage()
