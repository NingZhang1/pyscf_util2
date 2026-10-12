from pyscf import gto, scf, lib, symm
import numpy
from functools import reduce
from pyscf import lib
from pyscf import tools
import pyscf
from pyscf.symm import sph


def get_ao2L_spinor(mol, is_alpha=True):

    # n = mol.nao_2c()
    # nao = n // 2
    # JzLabel = numpy.diag(_get_Jz_AO(mol)[:n,:n])
    # res = numpy.zeros((nao, n), dtype=numpy.complex128)
    # loc_n = 0
    # loc_nao = 0
    # for i in range(mol.nbas):
    #     l = mol.bas_angular(i)
    #     nctr = mol.bas_nctr(i)
    #     if is_alpha:
    #         tmp = gto.sph2spinor_l(l)[0]
    #     else:
    #         tmp = gto.sph2spinor_l(l)[1]
    #     assert tmp.shape == (2 * l + 1, 4 * l + 2)
    #     res[loc_n : loc_n + 2 * l + 1, loc_nao : loc_nao + 4 * l + 2] = tmp
    #     loc_n += 2 * l + 1
    #     loc_nao += 4 * l + 2
    # return res

    ao2Lalpha_bench, ao2Lbeta_bench = sph.sph2spinor_coeff(mol)

    if is_alpha:
        return ao2Lalpha_bench
    else:
        return ao2Lbeta_bench


if __name__ == "__main__":

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
        verbose=5,
        charge=0,
        spin=0,
        symmetry="d2h",
    )

    # 获取每个 shell 的角动量、收缩数、原函数数
    shell_l = []
    shell_nctr = []
    shell_nprim = []
    for i in range(mol.nbas):
        l = mol.bas_angular(i)
        nctr = mol.bas_nctr(i)
        nprim = mol.bas_nprim(i)
        shell_l.append(l)
        shell_nctr.append(nctr)
        shell_nprim.append(nprim)

    # 每个 shell 产生的基函数个数（球谐，每个收缩 2l+1 个）
    shell_nfun = [(2 * l + 1) * nctr for l, nctr in zip(shell_l, shell_nctr)]

    # 计算每个 shell 的起始基函数索引
    start_idx = numpy.cumsum([0] + shell_nfun)[:-1]

    print("\n各壳层信息：")
    print("Shell  Atom  l  Symbol  nprim  nctr  nfun  Start")
    for i, (l, nprim, nctr, nfun, start) in enumerate(
        zip(shell_l, shell_nprim, shell_nctr, shell_nfun, start_idx)
    ):
        atom = mol.bas_atom(i)
        symbol = "spdfgh"[l] if l < len("spdfgh") else f"l={l}"
        print(
            f"{i:4d}  {atom:4d}  {l:2d}  {symbol:6s}  {nprim:5d}  {nctr:4d}  {nfun:4d}  {start:5d}"
        )

    # 打印每个标量基函数的角动量序列
    orb_ang = []
    for l, nfun in zip(shell_l, shell_nfun):
        orb_ang.extend([l] * nfun)

    print("\n每个标量基函数的角动量 (l) 序列：")
    print(orb_ang)

    # print(get_ao2L_spinor(mol))

    ao2Lalpha = get_ao2L_spinor(mol)
    ao2Lbeta = get_ao2L_spinor(mol, False)

    ovlp_scalar = mol.intor("int1e_ovlp")
    # print(ovlp_scalar)
    ovlp_spinor = (
        ao2Lalpha.conj().T @ ovlp_scalar @ ao2Lalpha
        + ao2Lbeta.conj().T @ ovlp_scalar @ ao2Lbeta
    )
    ovlp_spinor_bench = mol.intor("int1e_ovlp_spinor")

    print(numpy.allclose(ovlp_spinor, mol.intor("int1e_ovlp_spinor")))
