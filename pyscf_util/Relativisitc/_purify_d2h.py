from pyscf import gto, scf, lib, symm
import numpy
from functools import reduce
from pyscf import lib
from pyscf import tools
import pyscf


def label_orb_symm(mol, irrep_name, symm_orb, mo, s=None):
    nmo = mo.shape[1]
    if s is None:
        s = mol.intor_symmetric("int1e_ovlp")
    s_mo = numpy.dot(s, mo)
    norm = numpy.zeros((len(irrep_name), nmo))
    for i, csym in enumerate(symm_orb):
        moso = numpy.dot(csym.T.conj(), s_mo)
        ovlpso = reduce(numpy.dot, (csym.T.conj(), s, csym))
        try:
            s_moso = lib.cho_solve(ovlpso, moso)
        except numpy.linalg.LinAlgError:
            print("numpy.linalg.LinAlgError")
            ovlpso[numpy.diag_indices(csym.shape[1])] += 1e-12
            s_moso = lib.cho_solve(ovlpso, moso)
        norm[i] = numpy.einsum("ki,ki->i", moso.conj(), s_moso).real
    # print(norm)
    norm /= numpy.sum(norm, axis=0)  # for orbitals which are not normalized
    # iridx = numpy.argmax(norm, axis=0)
    # orbsym = numpy.asarray([irrep_name[i] for i in iridx])
    return norm


def purify_pes_d2h(mol, mf):

    # generate ao2L component #

    ao2Lalpha = get_ao2L_spinor(mol)
    ao2Lbeta = get_ao2L_spinor(mol, False)

    # generate mo_Lalpha, mo_Lbeta #

    n2c = mol.nao_2c()
    mo_Lalpha_coeff = ao2Lalpha @ mf.mo_coeff[:n2c, n2c:]
    mo_Lbeta_coeff = ao2Lbeta @ mf.mo_coeff[:n2c, n2c:]
    mo_Lalpha_coeff_real = mo_Lalpha_coeff.real
    mo_Lalpha_coeff_imag = mo_Lalpha_coeff.imag
    mo_Lbeta_coeff_real = mo_Lbeta_coeff.real
    mo_Lbeta_coeff_imag = mo_Lbeta_coeff.imag

    res = mf.mo_coeff[:, n2c:].copy()

    ovlp = mol.intor("int1e_ovlp")

    for i in range(0, n2c):
        norm_alpha_real = (
            mo_Lalpha_coeff_real[:, i].T @ ovlp @ mo_Lalpha_coeff_real[:, i]
        )
        norm_alpha_imag = (
            mo_Lalpha_coeff_imag[:, i].T @ ovlp @ mo_Lalpha_coeff_imag[:, i]
        )
        norm_beta_real = mo_Lbeta_coeff_real[:, i].T @ ovlp @ mo_Lbeta_coeff_real[:, i]
        norm_beta_imag = mo_Lbeta_coeff_imag[:, i].T @ ovlp @ mo_Lbeta_coeff_imag[:, i]

        norm_lst = [norm_alpha_real, norm_alpha_imag, norm_beta_real, norm_beta_imag]
        mo_coeff_lst = numpy.array(
            [
                mo_Lalpha_coeff_real[:, i],
                mo_Lalpha_coeff_imag[:, i],
                mo_Lbeta_coeff_real[:, i],
                mo_Lbeta_coeff_imag[:, i],
            ]
        )

        argmax_loc = numpy.argmax(numpy.array(norm_lst))

        # print(norm_lst)
        # print(argmax_loc)
        # assert argmax_loc in [0, 2]

        symm_loc = label_orb_symm(
            mol,
            mol.irrep_name,
            mol.symm_orb,
            mo_coeff_lst[argmax_loc].reshape(-1, 1),
            ovlp,
            # tol=1e-7/norm_lst[argmax_loc],
        )

        print(norm_lst[argmax_loc])
        print(symm_loc)
        # print(argmax_loc, mol.irrep_name[argmax_loc], norm_lst[argmax_loc])

        # print(
        #     f"MO {i} is assigned to irrep {mol.irrep_name[symm_loc.argmax()]} with norm {norm_lst[argmax_loc]}"
        # )



    # final check #


def purify_pes_c2v(mol, mf):

    # generate ao2L component #

    # #

    return


# purify for d2h and its subgroup #

if __name__ == "__main__":

    from Lalpha_ao2spinor import get_ao2L_spinor

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
        # basis="3-21g",
        basis="cc-pvdz-dk",
        verbose=5,
        charge=0,
        spin=0,
        symmetry="d2h",
    )

    print(get_ao2L_spinor(mol))

    ao2Lalpha = get_ao2L_spinor(mol)
    ao2Lbeta = get_ao2L_spinor(mol, False)

    # do scf #

    mol.build()
    mf = scf.dhf.RDHF(mol)
    mf.conv_tol = 1e-12
    mf.kernel()

    purify_pes_d2h(mol, mf)
