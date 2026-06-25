from pyscf import lib
import numpy
from functools import reduce
from pyscf import tools


from pyscf_util.Relativisitc._4C_Ints_dump import (
    _dump_2e,
    _dump_2e_outcore,
    _r_outcore_Coulomb,
    _r_outcore_Breit,
    _r_outcore_Gaunt,
)
from pyscf_util.Relativisitc._4C_Ints_dump import (
    coulomb_LLLL,
    coulomb_LLSS,
    coulomb_SSLL,
    coulomb_SSSS,
    breit_LSLS,
    breit_SLSL,
    breit_LSSL,
    breit_SLLS,
)


def FCIDUMP_Rela4C(
    mol,
    my_RDHF,
    with_breit=None,
    npes=None,
    filename="fcidump",
    mode="incore",
    orbsym_ID=None,
    IsComplex=True,
    tol=1e-8,
    debug=False,
    no_2e=False,
    retain_paldus_123=True,
    retain_paldus_456=True,
    only_dg=False,
):
    """Dump the relativistic 4-component integrals in FCIDUMP format

    Args:
        mol: a molecule object
        my_RDHF: a pyscf Restricted Dirac HF object
        filename: the filename of the FCIDUMP file

    Kwargs:
        with_breit: whether to include Breit term
        mode: the mode to dump the integrals
        debug: whether to return the integrals

    Returns:

    """

    assert mode in ["original", "incore", "outcore"]

    PREFIX = "RELA_4C_%d" % (numpy.random.randint(1, 19951201 + 1))

    n2c = mol.nao_2c()
    mo_coeff = my_RDHF.mo_coeff
    mo_coeff_mat = numpy.matrix(mo_coeff)

    if npes is None:
        npes = n2c

    mo_coeff_pes = mo_coeff_mat[:, n2c : n2c + npes]

    hcore = my_RDHF.get_hcore()
    h1e = reduce(numpy.dot, (mo_coeff_pes.H, hcore, mo_coeff_pes))

    n4c = 2 * n2c

    if with_breit is None:
        with_breit = my_RDHF.with_breit
        with_gaunt = my_RDHF.with_gaunt

    with_gaunt = my_RDHF.with_gaunt

    if with_breit:
        INT_LSLS_name = "int2e_breit_ssp1ssp2_spinor"
        INT_SLSL_name = "int2e_breit_sps1sps2_spinor"
        INT_LSSL_name = "int2e_breit_ssp1sps2_spinor"
        INT_SLLS_name = "int2e_breit_sps1ssp2_spinor"
    else:
        if with_gaunt:
            INT_LSLS_name = "int2e_ssp1ssp2_spinor"
            INT_SLSL_name = "int2e_sps1sps2_spinor"
            INT_LSSL_name = "int2e_ssp1sps2_spinor"
            INT_SLLS_name = "int2e_sps1ssp2_spinor"

    if mode == "original":
        int2e_res = numpy.zeros((n4c, n4c, n4c, n4c), dtype=numpy.complex128)
        c1 = 0.5 / lib.param.LIGHT_SPEED
        int2e_res[:n2c, :n2c, :n2c, :n2c] = mol.intor("int2e_spinor")  # LL LL
        tmp = mol.intor("int2e_spsp1_spinor") * c1**2
        int2e_res[n2c:, n2c:, :n2c, :n2c] = tmp  # SS LL
        int2e_res[:n2c, :n2c, n2c:, n2c:] = tmp.transpose(2, 3, 0, 1)  # LL SS
        int2e_res[n2c:, n2c:, n2c:, n2c:] = (
            mol.intor("int2e_spsp1spsp2_spinor") * c1**4
        )  # SS SS
        int2e_coulomb = lib.einsum("ijkl,ip->pjkl", int2e_res, mo_coeff_pes.conj())
        int2e_coulomb = lib.einsum("pjkl,jq->pqkl", int2e_coulomb, mo_coeff_pes)
        int2e_coulomb = lib.einsum("pqkl,kr->pqrl", int2e_coulomb, mo_coeff_pes.conj())
        int2e_coulomb = lib.einsum("pqrl,ls->pqrs", int2e_coulomb, mo_coeff_pes)
        if with_breit or with_gaunt:
            int2e_breit = numpy.zeros((n4c, n4c, n4c, n4c), dtype=numpy.complex128)
            ##### (LS|LS) and (SL|SL) #####
            tmp = mol.intor(INT_LSLS_name) * c1**2
            int2e_breit[:n2c, n2c:, :n2c, n2c:] = tmp
            tmp = mol.intor(INT_SLSL_name) * c1**2
            int2e_breit[n2c:, :n2c, n2c:, :n2c] = tmp
            ##### (LS|SL) and (SL|LS) #####
            tmp2 = mol.intor(INT_LSSL_name) * c1**2
            int2e_breit[:n2c, n2c:, n2c:, :n2c] = tmp2  # (LS|SL)
            tmp2 = mol.intor(INT_SLLS_name) * c1**2
            int2e_breit[n2c:, :n2c, :n2c, n2c:] = tmp2  # (SL|LS)
            ###############################
            int2e_breit = lib.einsum("ijkl,ip->pjkl", int2e_breit, mo_coeff_pes.conj())
            int2e_breit = lib.einsum("pjkl,jq->pqkl", int2e_breit, mo_coeff_pes)
            int2e_breit = lib.einsum("pqkl,kr->pqrl", int2e_breit, mo_coeff_pes.conj())
            int2e_breit = lib.einsum("pqrl,ls->pqrs", int2e_breit, mo_coeff_pes)
        else:
            int2e_breit = None
    elif mode == "incore":
        c1 = 0.5 / lib.param.LIGHT_SPEED
        ### LLLL part ###
        int2e_tmp = mol.intor("int2e_spinor")
        mo_coeff_L = mo_coeff_pes[:n2c, :]
        int2e_res = lib.einsum(
            "pqrs,pi,qj,rk,sl->ijkl",
            int2e_tmp,
            mo_coeff_L.conj(),
            mo_coeff_L,
            mo_coeff_L.conj(),
            mo_coeff_L,
        )
        ### SSLL part ###
        int2e_tmp = mol.intor("int2e_spsp1_spinor") * c1**2
        mo_coeff_S = mo_coeff_pes[n2c:, :]
        int2e_tmp = lib.einsum(
            "pqrs,pi,qj,rk,sl->ijkl",
            int2e_tmp,
            mo_coeff_S.conj(),
            mo_coeff_S,
            mo_coeff_L.conj(),
            mo_coeff_L,
        )
        int2e_res += int2e_tmp
        int2e_res += int2e_tmp.transpose(2, 3, 0, 1)
        ### SSSS part ###
        int2e_tmp = mol.intor("int2e_spsp1spsp2_spinor") * c1**4
        int2e_res += lib.einsum(
            "pqrs,pi,qj,rk,sl->ijkl",
            int2e_tmp,
            mo_coeff_S.conj(),
            mo_coeff_S,
            mo_coeff_S.conj(),
            mo_coeff_S,
        )
        int2e_coulomb = int2e_res
        if with_breit or with_gaunt:
            ### LSLS part ###
            int2e_tmp = mol.intor(INT_LSLS_name) * c1**2
            int2e_breit = lib.einsum(
                "pqrs,pi,qj,rk,sl->ijkl",
                int2e_tmp,
                mo_coeff_L.conj(),
                mo_coeff_S,
                mo_coeff_L.conj(),
                mo_coeff_S,
            )
            ### SLSL part ###
            int2e_tmp = int2e_tmp.conj().transpose(1, 0, 3, 2)
            int2e_breit += lib.einsum(
                "pqrs,pi,qj,rk,sl->ijkl",
                int2e_tmp,
                mo_coeff_S.conj(),
                mo_coeff_L,
                mo_coeff_S.conj(),
                mo_coeff_L,
            )
            ### LSSL part ###
            int2e_tmp = mol.intor(INT_LSSL_name) * c1**2
            int2e_breit += lib.einsum(
                "pqrs,pi,qj,rk,sl->ijkl",
                int2e_tmp,
                mo_coeff_L.conj(),
                mo_coeff_S,
                mo_coeff_S.conj(),
                mo_coeff_L,
            )
            ### SLLS part ###
            int2e_tmp = int2e_tmp.transpose(2, 3, 0, 1)
            int2e_breit += lib.einsum(
                "pqrs,pi,qj,rk,sl->ijkl",
                int2e_tmp,
                mo_coeff_S.conj(),
                mo_coeff_L,
                mo_coeff_L.conj(),
                mo_coeff_S,
            )
        else:
            int2e_breit = None
    elif mode == "outcore":

        # pyscf.ao2mo.r_outcore.general
        # raise NotImplementedError("outcore mode is not implemented yet")

        _r_outcore_Coulomb(mol, my_RDHF, npes, PREFIX)
        if with_breit:
            _r_outcore_Breit(mol, my_RDHF, npes, PREFIX)
        else:
            if with_gaunt:
                _r_outcore_Gaunt(mol, my_RDHF, npes, PREFIX)

        int2e_coulomb = None
        int2e_breit = None

    else:
        raise ValueError("Unknown mode %s" % mode)

    energy_core = mol.get_enuc()

    nmo = n2c // 2
    nelec = mol.nelectron
    ms = 0
    tol = 1e-8
    nuc = energy_core
    float_format = " %18.12E"

    if orbsym_ID is None:
        orbsym_ID = []
        for _ in range(nmo):
            orbsym_ID.append(0)
    else:
        orbsym_ID = orbsym_ID[nmo:]

    with open(filename, "w") as fout:  # 4-fold symmetry
        tools.fcidump.write_head(fout, nmo, nelec, ms, orbsym_ID)

        # output_format = float_format + float_format + ' %4d %4d %4d %4d\n'
        # if int2e_coulomb.ndim == 4:

        if not no_2e:
            if mode != "outcore":
                if debug:
                    _dump_2e(
                        fout,
                        int2e_coulomb,
                        int2e_breit,
                        with_breit,
                        IsComplex,
                        symmetry="s1",
                        tol=tol,
                        retain_paldus_123=retain_paldus_123,
                        retain_paldus_456=retain_paldus_456,
                        only_dg=only_dg,
                    )
                else:
                    _dump_2e(
                        fout,
                        int2e_coulomb,
                        int2e_breit,
                        with_breit,
                        IsComplex,
                        symmetry="s4",
                        tol=tol,
                        retain_paldus_123=retain_paldus_123,
                        retain_paldus_456=retain_paldus_456,
                        only_dg=only_dg,
                    )
            else:
                if debug:
                    _dump_2e_outcore(
                        fout,
                        npes,
                        PREFIX,
                        with_breit,
                        with_gaunt,
                        IsComplex,
                        symmetry="s1",
                        tol=tol,
                    )
                else:
                    _dump_2e_outcore(
                        fout,
                        npes,
                        PREFIX,
                        with_breit,
                        with_gaunt,
                        IsComplex,
                        symmetry="s4",
                        tol=tol,
                    )

        ############################################ DUMP E1 #############################################

        if IsComplex:
            output_format = float_format + float_format + " %4d %4d  0  0\n"
            for i in range(npes):
                # for j in range(n2c):
                for j in range(i + 1):
                    if only_dg and ((i // 2) != (j // 2)):
                        continue
                    if abs(h1e[i, j]) > tol:
                        fout.write(
                            output_format
                            % (h1e[i, j].real, h1e[i, j].imag, i + 1, j + 1)
                        )
            output_format = float_format + " 0.0  0  0  0  0\n"
            fout.write(output_format % nuc)
        else:
            output_format = float_format + " %4d %4d  0  0\n"
            for i in range(npes):
                # for j in range(n2c):
                for j in range(i + 1):
                    if only_dg and ((i // 2) != (j // 2)):
                        continue
                    if abs(h1e[i, j]) > tol:
                        fout.write(output_format % (h1e[i, j].real, i + 1, j + 1))
            output_format = float_format + " 0  0  0  0\n"
            fout.write(output_format % nuc)

    # clean

    if mode == "outcore":

        import os

        os.remove(coulomb_LLLL % PREFIX)
        os.remove(coulomb_LLSS % PREFIX)
        os.remove(coulomb_SSLL % PREFIX)
        os.remove(coulomb_SSSS % PREFIX)

        if with_breit or with_gaunt:
            os.remove(breit_LSLS % PREFIX)
            os.remove(breit_SLSL % PREFIX)
            os.remove(breit_LSSL % PREFIX)
            os.remove(breit_SLLS % PREFIX)

    if debug:
        return int2e_coulomb, int2e_breit
    else:
        return None, None


def FCIDUMP_Rela4C_SU2(
    mol,
    my_RDHF,
    with_breit=False,
    npes=None,
    filename="fcidump",
    mode="incore",
    debug=False,
):

    from pyscf_util.Relativisitc.double_group import (
        _atom_Jz_adapted,
        atm_d2h_symmetry_adapt_mo_coeff,
        _atom_spinor_spatial_parity,
    )

    ######## adapt the molecular orbitals to the Jz symmetry and tr symmetry ########

    mo_coeff = _atom_Jz_adapted(mol, my_RDHF.mo_coeff, my_RDHF.mo_energy, debug=debug)
    mo_pes = atm_d2h_symmetry_adapt_mo_coeff(mol, mo_coeff, debug)
    n2c = mol.nao_2c()
    mo_coeff[:, n2c:] = mo_pes

    mo_parity = _atom_spinor_spatial_parity(mol, mo_coeff, debug=debug)

    print("mo_parity = ", mo_parity)
    mo_parity = [x for id, x in enumerate(mo_parity) if id % 2 == 0]
    print("mo_parity = ", mo_parity)

    my_RDHF.mo_coeff = mo_coeff

    if npes is None:
        npes = n2c

    coulomb, breit = FCIDUMP_Rela4C(
        mol,
        my_RDHF,
        with_breit=with_breit,
        npes=npes,
        filename=filename,
        mode=mode,
        orbsym_ID=mo_parity,
        IsComplex=False,
        debug=debug,
    )

    if debug:
        FCIDUMP_Rela4C(
            mol,
            my_RDHF,
            with_breit=with_breit,
            filename=filename + "_complex",
            mode=mode,
            orbsym_ID=None,
            IsComplex=True,
            debug=debug,
        )

    return coulomb, breit, mo_parity, mo_coeff


if __name__ == "__main__":
    from pyscf import gto, scf

    # mol = gto.M(atom='H 0 0 0; H 0 0 1; O 0 1 0', basis='sto-3g', verbose=5)
    mol = gto.M(
        atom="F 0 0 0",
        basis="cc-pvdz-dk",
        verbose=5,
        charge=-1,
        spin=0,
        symmetry="d2h",
    )
    mol.build()
    mf = scf.dhf.RDHF(mol)
    mf.conv_tol = 1e-12
    mf.kernel()

    mf.with_breit = True
    mf.kernel()

    int2e1, breit_1 = FCIDUMP_Rela4C(
        mol,
        mf,
        with_breit=True,
        filename="FCIDUMP_4C_Breit",
        mode="original",
        debug=True,
    )

    int2e2, breit_2 = FCIDUMP_Rela4C(
        mol,
        mf,
        with_breit=True,
        filename="FCIDUMP_4C_incore",
        mode="incore",
        debug=True,
    )

    for i in range(mol.nao_2c(), mol.nao_2c() * 2):
        for j in range(mf.mo_coeff.shape[0]):
            if abs(mf.mo_coeff[j, i]) > 1e-6:
                print(
                    "%4d %4d %15.8f %15.8f"
                    % (
                        i - mol.nao_2c() + 1,
                        j,
                        mf.mo_coeff[j, i].real,
                        mf.mo_coeff[j, i].imag,
                    )
                )

    print("diff = ", numpy.linalg.norm(int2e1 - int2e2))

    if breit_1 is not None:
        print("breit diff = ", numpy.linalg.norm(breit_1 - breit_2))

    #### check breit term 4-fold symmetry ####

    nao = mol.nao

    for i in range(nao * 2):
        for j in range(nao * 2):
            for k in range(nao * 2):
                for l in range(nao * 2):
                    # print(breit_1[i,j,k,l], breit_1[j,i,l,k])
                    t1 = abs(breit_1[i, j, k, l] - breit_1[j, i, l, k].conj())
                    t2 = abs(breit_1[i, j, k, l] - breit_1[k, l, i, j].conj())
                    if t1 > 1e-8:
                        print("Breit 4-fold symmetry is not satisfied")
                        print(breit_1[i, j, k, l], breit_1[j, i, l, k])

    for i in range(nao):
        for j in range(nao):
            for k in range(nao):
                for l in range(nao):

                    t1 = abs(
                        int2e1[2 * i, 2 * j, 2 * k, 2 * l]
                        - int2e1[2 * i, 2 * j, 2 * l + 1, 2 * k + 1]
                    )
                    t2 = abs(
                        int2e1[2 * j + 1, 2 * i + 1, 2 * k, 2 * l]
                        - int2e1[2 * j + 1, 2 * i + 1, 2 * l + 1, 2 * k + 1]
                    )
                    t3 = abs(
                        int2e1[2 * i, 2 * j, 2 * k, 2 * l]
                        - int2e1[2 * j + 1, 2 * i + 1, 2 * k, 2 * l]
                    )
                    if t1 > 1e-8 or t2 > 1e-8 or t3 > 1e-8:
                        print("Coulomb AAAA group is not time-reversal symmetric")
                        print(
                            int2e1[2 * i, 2 * j, 2 * k, 2 * l],
                            int2e1[2 * i, 2 * j, 2 * l + 1, 2 * k + 1],
                            int2e1[2 * j + 1, 2 * i + 1, 2 * k, 2 * l],
                            int2e1[2 * j + 1, 2 * i + 1, 2 * l + 1, 2 * k + 1],
                        )

                    t1 = abs(
                        breit_1[2 * i, 2 * j, 2 * k, 2 * l]
                        + breit_1[2 * i, 2 * j, 2 * l + 1, 2 * k + 1]
                    )
                    t2 = abs(
                        breit_1[2 * j + 1, 2 * i + 1, 2 * k, 2 * l]
                        + breit_1[2 * j + 1, 2 * i + 1, 2 * l + 1, 2 * k + 1]
                    )
                    t3 = abs(
                        breit_1[2 * i, 2 * j, 2 * k, 2 * l]
                        + breit_1[2 * j + 1, 2 * i + 1, 2 * k, 2 * l]
                    )
                    if t1 > 1e-8 or t2 > 1e-8 or t3 > 1e-8:
                        print("Breit AAAA group is not time-reversal symmetric")
                        print(
                            breit_1[2 * i, 2 * j, 2 * k, 2 * l + 1],
                            -breit_1[2 * i, 2 * j, 2 * l + 1, 2 * k + 1],
                            -breit_1[2 * j + 1, 2 * i + 1, 2 * k, 2 * l],
                            breit_1[2 * j + 1, 2 * i + 1, 2 * l + 1, 2 * k + 1],
                        )

                    t1 = abs(
                        int2e1[2 * i, 2 * j, 2 * k, 2 * l + 1]
                        + int2e1[2 * i, 2 * j, 2 * l, 2 * k + 1]
                    )
                    t2 = abs(
                        int2e1[2 * j + 1, 2 * i + 1, 2 * k, 2 * l + 1]
                        + int2e1[2 * j + 1, 2 * i + 1, 2 * l, 2 * k + 1]
                    )
                    t3 = abs(
                        int2e1[2 * i, 2 * j, 2 * k, 2 * l + 1]
                        - int2e1[2 * j + 1, 2 * i + 1, 2 * k, 2 * l + 1]
                    )
                    if t1 > 1e-8 or t2 > 1e-8 or t3 > 1e-8:
                        print("Coulomb AAAB group is not time-reversal symmetric")
                        print(
                            int2e1[2 * i, 2 * j, 2 * k, 2 * l + 1],
                            -int2e1[2 * i, 2 * j, 2 * l, 2 * k + 1],
                            int2e1[2 * j + 1, 2 * i + 1, 2 * k, 2 * l + 1],
                            -int2e1[2 * j + 1, 2 * i + 1, 2 * l, 2 * k + 1],
                        )

                    t1 = abs(
                        breit_1[2 * i, 2 * j, 2 * k, 2 * l + 1]
                        - breit_1[2 * i, 2 * j, 2 * l, 2 * k + 1]
                    )
                    t2 = abs(
                        breit_1[2 * j + 1, 2 * i + 1, 2 * k, 2 * l + 1]
                        - breit_1[2 * j + 1, 2 * i + 1, 2 * l, 2 * k + 1]
                    )
                    t3 = abs(
                        breit_1[2 * i, 2 * j, 2 * k, 2 * l + 1]
                        + breit_1[2 * j + 1, 2 * i + 1, 2 * k, 2 * l + 1]
                    )
                    if t1 > 1e-8 or t2 > 1e-8 or t3 > 1e-8:
                        print("Breit AAAB group is not time-reversal symmetric")
                        print(
                            breit_1[2 * i, 2 * j, 2 * k, 2 * l + 1],
                            breit_1[2 * i, 2 * j, 2 * l, 2 * k + 1],
                            -breit_1[2 * j + 1, 2 * i + 1, 2 * k, 2 * l + 1],
                            -breit_1[2 * j + 1, 2 * i + 1, 2 * l, 2 * k + 1],
                        )

                    t1 = abs(
                        int2e1[2 * i, 2 * j + 1, 2 * k, 2 * l + 1]
                        + int2e1[2 * i, 2 * j + 1, 2 * l, 2 * k + 1]
                    )
                    t2 = abs(
                        int2e1[2 * j, 2 * i + 1, 2 * k, 2 * l + 1]
                        + int2e1[2 * j, 2 * i + 1, 2 * l, 2 * k + 1]
                    )
                    t3 = abs(
                        int2e1[2 * i, 2 * j + 1, 2 * k, 2 * l + 1]
                        + int2e1[2 * j, 2 * i + 1, 2 * k, 2 * l + 1]
                    )
                    if t1 > 1e-8 or t2 > 1e-8 or t3 > 1e-8:
                        print("Coulomb ABAB group is not time-reversal symmetric")
                        print(
                            int2e1[2 * i, 2 * j + 1, 2 * k, 2 * l + 1],
                            -int2e1[2 * i, 2 * j + 1, 2 * l, 2 * k + 1],
                            -int2e1[2 * j, 2 * i + 1, 2 * k, 2 * l + 1],
                            int2e1[2 * j, 2 * i + 1, 2 * l, 2 * k + 1],
                        )

                    t1 = abs(
                        breit_1[2 * i, 2 * j + 1, 2 * k, 2 * l + 1]
                        - breit_1[2 * i, 2 * j + 1, 2 * l, 2 * k + 1]
                    )
                    t2 = abs(
                        breit_1[2 * j, 2 * i + 1, 2 * k, 2 * l + 1]
                        - breit_1[2 * j, 2 * i + 1, 2 * l, 2 * k + 1]
                    )
                    t3 = abs(
                        breit_1[2 * i, 2 * j + 1, 2 * k, 2 * l + 1]
                        - breit_1[2 * j, 2 * i + 1, 2 * k, 2 * l + 1]
                    )
                    if t1 > 1e-8 or t2 > 1e-8 or t3 > 1e-8:
                        print("Breit ABAB group is not time-reversal symmetric")
                        print(
                            breit_1[2 * i, 2 * j + 1, 2 * k, 2 * l + 1],
                            breit_1[2 * i, 2 * j + 1, 2 * l, 2 * k + 1],
                            breit_1[2 * j, 2 * i + 1, 2 * k, 2 * l + 1],
                            breit_1[2 * j, 2 * i + 1, 2 * l, 2 * k + 1],
                        )

                    t1 = abs(
                        int2e1[2 * i + 1, 2 * j, 2 * k, 2 * l + 1]
                        + int2e1[2 * i + 1, 2 * j, 2 * l, 2 * k + 1]
                    )
                    t2 = abs(
                        int2e1[2 * j + 1, 2 * i, 2 * k, 2 * l + 1]
                        + int2e1[2 * j + 1, 2 * i, 2 * l, 2 * k + 1]
                    )
                    t3 = abs(
                        int2e1[2 * i + 1, 2 * j, 2 * k, 2 * l + 1]
                        + int2e1[2 * j + 1, 2 * i, 2 * k, 2 * l + 1]
                    )
                    if t1 > 1e-8 or t2 > 1e-8 or t3 > 1e-8:
                        print("Coulomb BAAB group is not time-reversal symmetric")
                        print(
                            int2e1[2 * i + 1, 2 * j, 2 * k, 2 * l + 1],
                            -int2e1[2 * i + 1, 2 * j, 2 * l, 2 * k + 1],
                            -int2e1[2 * j + 1, 2 * i, 2 * k, 2 * l + 1],
                            int2e1[2 * j + 1, 2 * i, 2 * l, 2 * k + 1],
                        )

                    t1 = abs(
                        breit_1[2 * i + 1, 2 * j, 2 * k, 2 * l + 1]
                        - breit_1[2 * i + 1, 2 * j, 2 * l, 2 * k + 1]
                    )
                    t2 = abs(
                        breit_1[2 * j + 1, 2 * i, 2 * k, 2 * l + 1]
                        - breit_1[2 * j + 1, 2 * i, 2 * l, 2 * k + 1]
                    )
                    t3 = abs(
                        breit_1[2 * i + 1, 2 * j, 2 * k, 2 * l + 1]
                        - breit_1[2 * j + 1, 2 * i, 2 * k, 2 * l + 1]
                    )
                    if t1 > 1e-8 or t2 > 1e-8 or t3 > 1e-8:
                        print("Breit BAAB group is not time-reversal symmetric")
                        print(
                            breit_1[2 * i + 1, 2 * j, 2 * k, 2 * l + 1],
                            breit_1[2 * i + 1, 2 * j, 2 * l, 2 * k + 1],
                            breit_1[2 * j + 1, 2 * i, 2 * k, 2 * l + 1],
                            breit_1[2 * j + 1, 2 * i, 2 * l, 2 * k + 1],
                        )

    ##### check out_core mode #####

    def view(h5file, dataname="eri_mo"):
        import h5py

        f5 = h5py.File(h5file, "r")
        print("dataset %s, shape %s" % (str(f5.keys()), str(f5[dataname].shape)))
        f5.close()

    n2c = 2 * nao

    _r_outcore_Coulomb(mol, mf, n2c, prefix="F", max_memory=5, ioblk_size=2)
    _r_outcore_Breit(mol, mf, n2c, prefix="F", max_memory=5, ioblk_size=2)

    view(coulomb_LLLL % "F")
    view(coulomb_LLSS % "F")
    view(coulomb_SSLL % "F")
    view(coulomb_SSSS % "F")

    view(breit_LSLS % "F")
    view(breit_SLSL % "F")
    view(breit_LSSL % "F")
    view(breit_SLLS % "F")

    ### check whether it is correct ###

    import h5py

    n = mol.nao_2c()
    c1 = 0.5 / lib.param.LIGHT_SPEED

    feri_coulomb_LLLL = h5py.File(coulomb_LLLL % "F", "r")
    feri_coulomb_LLSS = h5py.File(coulomb_LLSS % "F", "r")
    feri_coulomb_SSLL = h5py.File(coulomb_SSLL % "F", "r")
    feri_coulomb_SSSS = h5py.File(coulomb_SSSS % "F", "r")

    eri_coulomb = numpy.array(feri_coulomb_LLLL["eri_mo"]).reshape(n, n, n, n)
    eri_coulomb += numpy.array(feri_coulomb_LLSS["eri_mo"]).reshape(n, n, n, n) * c1**2
    eri_coulomb += numpy.array(feri_coulomb_SSLL["eri_mo"]).reshape(n, n, n, n) * c1**2
    eri_coulomb += numpy.array(feri_coulomb_SSSS["eri_mo"]).reshape(n, n, n, n) * c1**4

    feri_breit_LSLS = h5py.File(breit_LSLS % "F", "r")
    feri_breit_SLSL = h5py.File(breit_SLSL % "F", "r")
    feri_breit_LSSL = h5py.File(breit_LSSL % "F", "r")
    feri_breit_SLLS = h5py.File(breit_SLLS % "F", "r")

    eri_breit = numpy.array(feri_breit_LSLS["eri_mo"]).reshape(n, n, n, n)
    eri_breit += numpy.array(feri_breit_SLSL["eri_mo"]).reshape(n, n, n, n)
    eri_breit += numpy.array(feri_breit_LSSL["eri_mo"]).reshape(n, n, n, n)
    eri_breit += numpy.array(feri_breit_SLLS["eri_mo"]).reshape(n, n, n, n)
    eri_breit *= c1**2

    print("eri_coulomb diff = ", numpy.linalg.norm(eri_coulomb - int2e1))
    print("eri_breit diff = ", numpy.linalg.norm(eri_breit - breit_1))

    int2e2, breit_2 = FCIDUMP_Rela4C(
        mol,
        mf,
        with_breit=True,
        filename="FCIDUMP_4C_outcore_1",
        mode="outcore",
        debug=True,
    )
    int2e2, breit_2 = FCIDUMP_Rela4C(
        mol,
        mf,
        with_breit=True,
        filename="FCIDUMP_4C_outcore_2",
        mode="outcore",
        debug=False,
    )

    feri_coulomb_LLLL.close()
    feri_coulomb_LLSS.close()
    feri_coulomb_SSLL.close()
    feri_coulomb_SSSS.close()

    feri_breit_LSLS.close()
    feri_breit_SLSL.close()
    feri_breit_LSSL.close()
    feri_breit_SLLS.close()

    int_coulomb, int_breit, mo_parity, mo_coeff_adapted = FCIDUMP_Rela4C_SU2(
        mol, mf, with_breit=True, filename="fcidump_adapted", mode="incore", debug=True
    )
