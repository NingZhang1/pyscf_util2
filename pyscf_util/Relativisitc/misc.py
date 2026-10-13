"""Wigner rotation D-matrix for J=1/2"""

import numpy
from scipy.special import factorial


def dmatrix_half(beta):
    """Wigner small-d matrix (in z-y-z convention) for J=1/2"""
    c = numpy.cos(beta / 2)
    s = numpy.sin(beta / 2)
    mat = numpy.array(((c, -s), (s, c)), dtype=numpy.complex128)
    return mat


def Dmatrix_half(alpha, beta, gamma):
    """Wigner rotation D-matrix

    D_{mm'} = <lm|R(alpha,beta,gamma)|lm'>
    alpha, beta, gamma are Euler angles (in z-y-z convention)
    l = 1/2

    Kwargs:
        reorder_p (bool): Whether to put the p functions in the (x,y,z) order.
    """
    d = dmatrix_half(beta)
    ms = numpy.asarray([0.5, -0.5])  # |alpha>, |beta>
    D = numpy.einsum(
        "i,ij,j->ij", numpy.exp(-1j * alpha * ms), d, numpy.exp(-1j * gamma * ms)
    )

    return D


import numpy as np
from pyscf.symm import geom


def mol_principal_axes(mol):
    """
    返回分子的主惯性轴。
    输出：(I1, I2, I3), (e1, e2, e3)
    其中 e_i 为对应的单位主轴向量（列向量形式），转动惯量按升序排列。
    """
    coords = mol.atom_coords(unit="Bohr")  # shape (natm, 3)
    masses = np.asarray(mol.atom_mass_list())  # 原子质量
    # 质心
    com = (masses[:, None] * coords).sum(axis=0) / masses.sum()
    r = coords - com
    # 惯性张量
    I = np.zeros((3, 3))
    for m, (x, y, z) in zip(masses, r):
        I += m * np.array(
            [
                [y * y + z * z, -x * y, -x * z],
                [-x * y, x * x + z * z, -y * z],
                [-x * z, -y * z, x * x + y * y],
            ]
        )
    # 特征分解（对称矩阵，eigh 保证正交矩阵）
    vals, vecs = np.linalg.eigh(I)  # vals 升序
    return vals, vecs  # vecs 的列向量即三个主轴


def handedness(mol):
    """
    判断分子对称性主轴是否构成右手系。
    返回 (是否右手系, det, 转动惯量)
    """
    vals, vecs = mol_principal_axes(mol)
    _, _, vecs = geom.detect_symm(mol._atom)
    det = float(np.linalg.det(vecs))
    # 若转动惯量有简并（如线性分子、球陀螺），主轴不唯一，结果没有绝对意义
    degenerate = np.any(np.diff(vals) < 1e-8)
    return det > 0, det, vals, degenerate


if __name__ == "__main__":

    from pyscf.symm.Dmatrix import Dmatrix
    from pyscf import gto

    c1 = [1, 0, 0]
    c2 = [-1, 0, 0]

    Ebar = [0, 0, 2 * numpy.pi]
    RX = [numpy.pi, numpy.pi, 0]
    RY = [0, numpy.pi, 0]
    RZ = [0, 0, numpy.pi]

    print(Dmatrix_half(Ebar[0], Ebar[1], Ebar[2]))
    print(Dmatrix_half(RX[0], RX[1], RX[2]))
    print(Dmatrix_half(RY[0], RY[1], RY[2]))
    print(Dmatrix_half(RZ[0], RZ[1], RZ[2]))

    print(Dmatrix(1, RX[0], RX[1], RX[2], True))
    print(Dmatrix(1, RY[0], RY[1], RY[2], True))
    print(Dmatrix(1, RZ[0], RZ[1], RZ[2], True))

    print(numpy.diag(Dmatrix(2, RX[0], RX[1], RX[2], True)))
    print(numpy.diag(Dmatrix(2, RY[0], RY[1], RY[2], True)))
    print(numpy.diag(Dmatrix(2, RZ[0], RZ[1], RZ[2], True)))

    print(numpy.diag(Dmatrix(3, RX[0], RX[1], RX[2], True)))
    print(numpy.diag(Dmatrix(3, RY[0], RY[1], RY[2], True)))
    print(numpy.diag(Dmatrix(3, RZ[0], RZ[1], RZ[2], True)))

    print(numpy.diag(Dmatrix(4, RX[0], RX[1], RX[2], True)))
    print(numpy.diag(Dmatrix(4, RY[0], RY[1], RY[2], True)))
    print(numpy.diag(Dmatrix(4, RZ[0], RZ[1], RZ[2], True)))

    mol = gto.M(
        atom="""
C   0.000000   0.000000   0.669500
C   0.000000   0.000000  -0.669500
H   0.928900   0.000000   1.232100
H  -0.928900   0.000000   1.232100
H   0.928900   0.000000  -1.232100
H  -0.928900   0.000000  -1.232100
""",
        # basis="unc-cc-pvdz-dk",
        basis="3-21g",
        # basis="cc-pvdz-dk",
        verbose=5,
        charge=0,
        spin=0,
        symmetry="d2h",
        # unit='bohr'
    )

    print(mol_principal_axes(mol))

    right, det, vals, degen = handedness(mol)
    print("转动惯量:", vals)
    print("det =", det)
    print("右手系" if right else "左手系")
    print("简并:", degen)

    _, _, vecs = geom.detect_symm(mol._atom)

    print(vecs)
    print(float(np.linalg.det(vecs)))
