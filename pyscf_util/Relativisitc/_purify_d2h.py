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
    norm_square = numpy.zeros((len(irrep_name), nmo))
    mo_coeff_sym = []
    for i, csym in enumerate(symm_orb):
        # print("i, csym", i, csym)
        moso = numpy.dot(csym.T.conj(), s_mo)
        ovlpso = reduce(numpy.dot, (csym.T.conj(), s, csym))
        try:
            s_moso = lib.cho_solve(ovlpso, moso)
        except numpy.linalg.LinAlgError:
            print("numpy.linalg.LinAlgError")
            ovlpso[numpy.diag_indices(csym.shape[1])] += 1e-12
            s_moso = lib.cho_solve(ovlpso, moso)
        # print(csym.shape)
        # print(s_moso.shape)
        # exit(1)
        norm_square[i] = numpy.einsum(
            "ki,ki->i", moso.conj(), s_moso
        ).real  # norm_square_indeed
        mo_coeff_sym.append((csym @ s_moso).reshape(-1))
    # print(norm)
    norm_square /= numpy.sum(
        norm_square, axis=0
    )  # for orbitals which are not normalized
    # 归一化 mo_coeff_sym #
    mo_coeff_sym = numpy.asarray(mo_coeff_sym).T
    mo_coeff_sym_ovlp = mo_coeff_sym.T @ s @ mo_coeff_sym
    mo_coeff_sym_offdg = mo_coeff_sym_ovlp - numpy.diag(numpy.diag(mo_coeff_sym_ovlp))
    assert numpy.allclose(mo_coeff_sym_offdg, 0.0)
    norm_mo_coeff_sym = numpy.sqrt(numpy.diag(mo_coeff_sym_ovlp))
    for i in range(mo_coeff_sym.shape[1]):
        mo_coeff_sym[:, i] /= norm_mo_coeff_sym[i]
    return norm_square, mo_coeff_sym


def purify_pes_d2h(mol, mf, pure_tol=1e-10):

    IRREP_NAME_D2H_FULL = ["Ag", "B1g", "B2g", "B3g", "Au", "B1u", "B2u", "B3u"]

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

    for i in range(0, n2c, 2):  # loop over unbarred only #

        print("---------------------------------------")
        print("i = ", i)

        ##################################################################################
        ## Stage I 处理不需要混合，只需要交换就能解决的纯化问题 ##
        ## \{|\alpha\rangle, |\beta\rangle\} --> \{|\beta\rangle, -|\alpha\rangle\}  ##
        ##################################################################################

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

        ### judge whether it is purified ###

        argmax_loc = numpy.argmax(numpy.array(norm_lst))

        symm_comp_norm_square, _ = label_orb_symm(
            mol,
            mol.irrep_name,
            mol.symm_orb,
            mo_coeff_lst[argmax_loc].reshape(-1, 1),
            ovlp,
            # tol=1e-7/norm_lst[argmax_loc],
        )

        argmax_symm_comp = numpy.argmax(numpy.array(symm_comp_norm_square))
        argmax_symm_comp_value = symm_comp_norm_square[argmax_symm_comp]
        argmax_symm_comp_irrep_name = mol.irrep_name[argmax_symm_comp]

        print("argmax_loc ", argmax_loc)
        print(
            "symm comp %s %20.16f"
            % (argmax_symm_comp_irrep_name, argmax_symm_comp_value[0])
        )

        is_pure = True

        if argmax_symm_comp_value < (1 - pure_tol):
            is_pure = False

        # print("is_pure ", is_pure)

        need_update = False

        # if is_pure:
        if argmax_loc == 0:  # alpha real
            if argmax_symm_comp_irrep_name not in ["Ag", "Au"]:
                is_pure = False
                if argmax_symm_comp_irrep_name not in ["B1g", "B1u"]:
                    need_update = True
                    res[:, i], res[:, i + 1] = res[:, i + 1], res[:, i]
                    res[:, i + 1] *= -1.0
        elif argmax_loc == 1:  # alpha imag
            if argmax_symm_comp_irrep_name not in ["B1g", "B1u"]:
                is_pure = False
                if argmax_symm_comp_irrep_name not in ["Ag", "Au"]:
                    need_update = True
                    res[:, i], res[:, i + 1] = res[:, i + 1], res[:, i]
                    res[:, i + 1] *= -1.0
        elif argmax_loc == 2:  # beta real
            if argmax_symm_comp_irrep_name not in ["B2g", "B2u"]:
                is_pure = False
                if argmax_symm_comp_irrep_name not in ["B3g", "B3u"]:
                    need_update = True
                    res[:, i], res[:, i + 1] = res[:, i + 1], res[:, i]
                    res[:, i + 1] *= -1.0
        elif argmax_loc == 3:  # beta imag
            if argmax_symm_comp_irrep_name not in ["B3g", "B3u"]:
                is_pure = False
                if argmax_symm_comp_irrep_name not in ["B2g", "B2u"]:
                    need_update = True
                    res[:, i], res[:, i + 1] = res[:, i + 1], res[:, i]
                    res[:, i + 1] *= -1.0
        else:
            raise RuntimeError

        print("is_pure ", is_pure)

        if is_pure:
            continue

        # update #

        if need_update:
            mo_Lalpha_coeff = ao2Lalpha @ res[:n2c, :]
            mo_Lbeta_coeff = ao2Lbeta @ res[:n2c, :]
            mo_Lalpha_coeff_real = mo_Lalpha_coeff.real
            mo_Lalpha_coeff_imag = mo_Lalpha_coeff.imag
            mo_Lbeta_coeff_real = mo_Lbeta_coeff.real
            mo_Lbeta_coeff_imag = mo_Lbeta_coeff.imag

        ##################################################################################
        ## Stage II 只需要转动basis的相位就能解决的纯化问题 ##
        ##################################################################################

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

        ### judge whether it is purified ###

        argmax_loc = numpy.argmax(numpy.array(norm_lst))

        symm_comp_norm_square, symm_comp_mocoeff = label_orb_symm(
            mol,
            mol.irrep_name,
            mol.symm_orb,
            mo_coeff_lst[argmax_loc].reshape(-1, 1),
            ovlp,
            # tol=1e-7/norm_lst[argmax_loc],
        )

        argmax_symm_comp = numpy.argmax(numpy.array(symm_comp_norm_square))
        argmax_symm_comp_value = symm_comp_norm_square[argmax_symm_comp]
        argmax_symm_comp_irrep_name = mol.irrep_name[argmax_symm_comp]

        symm_comp_norm_square_full = {}
        symm_comp_mocoeff_full = {}
        for x, name in enumerate(mol.irrep_name):
            symm_comp_norm_square_full[name] = symm_comp_norm_square[x]
            symm_comp_mocoeff_full[name] = symm_comp_mocoeff[x]
        for x in IRREP_NAME_D2H_FULL:
            if x not in symm_comp_norm_square_full:
                symm_comp_norm_square_full[x] = 0.0
                symm_comp_mocoeff_full[name] = None

        print("argmax_loc ", argmax_loc)
        print(
            "symm comp %s %20.16f"
            % (argmax_symm_comp_irrep_name, argmax_symm_comp_value[0])
        )

        is_pure = True

        if argmax_symm_comp_value < (1 - pure_tol):
            is_pure = False

        need_update = False

        if is_pure:
            if argmax_loc == 0:  # alpha real
                if argmax_symm_comp_irrep_name not in ["Ag", "Au"]:
                    assert argmax_symm_comp_irrep_name in ["B1g", "B1u"]
                    res[:, i] *= 1.0j
                    res[:, i + 1] *= -1.0j
                    need_update = True
            elif argmax_loc == 1:  # alpha imag
                if argmax_symm_comp_irrep_name not in ["B1g", "B1u"]:
                    assert argmax_symm_comp_irrep_name in ["Ag", "Au"]
                    res[:, i] *= 1.0j
                    res[:, i + 1] *= -1.0j
                    need_update = True
            elif argmax_loc == 2:  # beta real
                if argmax_symm_comp_irrep_name not in ["B2g", "B2u"]:
                    assert argmax_symm_comp_irrep_name in ["B3g", "B3u"]
                    res[:, i] *= 1.0j
                    res[:, i + 1] *= -1.0j
                    need_update = True
            elif argmax_loc == 3:  # beta imag
                if argmax_symm_comp_irrep_name not in ["B3g", "B3u"]:
                    assert argmax_symm_comp_irrep_name in ["B2g", "B2u"]
                    res[:, i] *= 1.0j
                    res[:, i + 1] *= -1.0j
                    need_update = True
            else:
                raise RuntimeError

        else:
            print(symm_comp_norm_square_full)

            if argmax_loc in [0, 1]:  # alpha real/imag
                if argmax_symm_comp_irrep_name in ["Ag", "B1g"]:
                    sum_comp = (
                        symm_comp_norm_square_full["Ag"]
                        + symm_comp_norm_square_full["B1g"]
                    )
                elif argmax_symm_comp_irrep_name in ["Au", "B1u"]:
                    sum_comp = (
                        symm_comp_norm_square_full["Au"]
                        + symm_comp_norm_square_full["B1u"]
                    )
                else:
                    raise RuntimeError
            elif argmax_loc in [2, 3]:  # beta real/imag
                if argmax_symm_comp_irrep_name in ["B2g", "B3g"]:
                    sum_comp = (
                        symm_comp_norm_square_full["B2g"]
                        + symm_comp_norm_square_full["B3g"]
                    )
                elif argmax_symm_comp_irrep_name in ["B2u", "B3u"]:
                    sum_comp = (
                        symm_comp_norm_square_full["B2u"]
                        + symm_comp_norm_square_full["B3u"]
                    )
                else:
                    raise RuntimeError
            else:
                raise RuntimeError

            if sum_comp < (1 - pure_tol):
                is_pure = False
            else:
                is_pure = True
                need_update = True

                #######################################
                ## 求解相位 ##
                #######################################

                if argmax_loc == 0:  # alpha real
                    if argmax_symm_comp_irrep_name in ["Ag", "B1g"]:
                        c = numpy.sqrt(symm_comp_norm_square_full["Ag"])
                        d = -numpy.sqrt(symm_comp_norm_square_full["B1g"])
                    elif argmax_symm_comp_irrep_name in ["Au", "B1u"]:
                        c = numpy.sqrt(symm_comp_norm_square_full["Au"])
                        d = -numpy.sqrt(symm_comp_norm_square_full["B1u"])
                    else:
                        raise RuntimeError
                elif argmax_loc == 1:  # alpha imag
                    if argmax_symm_comp_irrep_name in ["Ag", "B1g"]:
                        c = numpy.sqrt(symm_comp_norm_square_full["B1g"])
                        d = numpy.sqrt(symm_comp_norm_square_full["Ag"])
                    elif argmax_symm_comp_irrep_name in ["Au", "B1u"]:
                        c = numpy.sqrt(symm_comp_norm_square_full["B1u"])
                        d = numpy.sqrt(symm_comp_norm_square_full["Au"])
                    else:
                        raise RuntimeError
                elif argmax_loc == 2:
                    if argmax_symm_comp_irrep_name in ["B2g", "B3g"]:
                        c = numpy.sqrt(symm_comp_norm_square_full["B2g"])
                        d = -numpy.sqrt(symm_comp_norm_square_full["B3g"])
                    elif argmax_symm_comp_irrep_name in ["B2u", "B3u"]:
                        c = numpy.sqrt(symm_comp_norm_square_full["B2u"])
                        d = -numpy.sqrt(symm_comp_norm_square_full["B3u"])
                    else:
                        raise RuntimeError
                elif argmax_loc == 3:
                    if argmax_symm_comp_irrep_name in ["B2g", "B3g"]:
                        c = numpy.sqrt(symm_comp_norm_square_full["B3g"])
                        d = numpy.sqrt(symm_comp_norm_square_full["B2g"])
                    elif argmax_symm_comp_irrep_name in ["B2u", "B3u"]:
                        c = numpy.sqrt(symm_comp_norm_square_full["B3u"])
                        d = numpy.sqrt(symm_comp_norm_square_full["B2u"])
                    else:
                        raise RuntimeError
                else:
                    raise RuntimeError

                res[:, i] *= c + d * (1.0j)
                res[:, i + 1] *= c - d * (1.0j)
                need_update = True

        # update #

        if need_update:
            mo_Lalpha_coeff = ao2Lalpha @ res[:n2c, :]
            mo_Lbeta_coeff = ao2Lbeta @ res[:n2c, :]
            mo_Lalpha_coeff_real = mo_Lalpha_coeff.real
            mo_Lalpha_coeff_imag = mo_Lalpha_coeff.imag
            mo_Lbeta_coeff_real = mo_Lbeta_coeff.real
            mo_Lbeta_coeff_imag = mo_Lbeta_coeff.imag

        print("is_pure ", is_pure)

        if is_pure:
            continue

        ##################################################################################
        ## Stage III 需要转动 alpha beta 才能解决的纯化问题 ##
        ##################################################################################

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

        ### judge whether it is purified ###

        argmax_loc = numpy.argmax(numpy.array(norm_lst))

        symm_comp_norm_square, symm_comp_mocoeff = label_orb_symm(
            mol,
            mol.irrep_name,
            mol.symm_orb,
            mo_coeff_lst[argmax_loc].reshape(-1, 1),
            ovlp,
            # tol=1e-7/norm_lst[argmax_loc],
        )

        argmax_symm_comp = numpy.argmax(numpy.array(symm_comp_norm_square))
        argmax_symm_comp_value = symm_comp_norm_square[argmax_symm_comp]
        argmax_symm_comp_irrep_name = mol.irrep_name[argmax_symm_comp]

        symm_comp_norm_square_full = {}
        symm_comp_mocoeff_full = {}
        for x, name in enumerate(mol.irrep_name):
            symm_comp_norm_square_full[name] = symm_comp_norm_square[x]
            symm_comp_mocoeff_full[name] = symm_comp_mocoeff[x]
        for x in IRREP_NAME_D2H_FULL:
            if x not in symm_comp_norm_square_full:
                symm_comp_norm_square_full[x] = 0.0
                symm_comp_mocoeff_full[name] = None

        print("argmax_loc ", argmax_loc)
        print(
            "symm comp %s %20.16f"
            % (argmax_symm_comp_irrep_name, argmax_symm_comp_value[0])
        )

        #######################################
        ## 求解相位 ##
        #######################################

        if argmax_loc == 0:  # alpha real
            if argmax_symm_comp_irrep_name in ["Ag", "B1g"]:
                c = numpy.sqrt(symm_comp_norm_square_full["Ag"])
                d = -numpy.sqrt(symm_comp_norm_square_full["B1g"])
                e = numpy.sqrt(symm_comp_norm_square_full["B2g"])
                f = -numpy.sqrt(symm_comp_norm_square_full["B3g"])
            elif argmax_symm_comp_irrep_name in ["Au", "B1u"]:
                c = numpy.sqrt(symm_comp_norm_square_full["Au"])
                d = -numpy.sqrt(symm_comp_norm_square_full["B1u"])
                e = numpy.sqrt(symm_comp_norm_square_full["B2u"])
                f = -numpy.sqrt(symm_comp_norm_square_full["B3u"])
            else:
                raise RuntimeError
        elif argmax_loc == 1:  # alpha imag
            if argmax_symm_comp_irrep_name in ["Ag", "B1g"]:
                c = numpy.sqrt(symm_comp_norm_square_full["B1g"])
                d = numpy.sqrt(symm_comp_norm_square_full["Ag"])
                e = numpy.sqrt(symm_comp_norm_square_full["B3g"])
                f = numpy.sqrt(symm_comp_norm_square_full["B2g"])
            elif argmax_symm_comp_irrep_name in ["Au", "B1u"]:
                c = numpy.sqrt(symm_comp_norm_square_full["B1u"])
                d = numpy.sqrt(symm_comp_norm_square_full["Au"])
                e = numpy.sqrt(symm_comp_norm_square_full["B3u"])
                f = numpy.sqrt(symm_comp_norm_square_full["B2u"])
            else:
                raise RuntimeError
        elif argmax_loc == 2:
            if argmax_symm_comp_irrep_name in ["B2g", "B3g"]:
                c = numpy.sqrt(symm_comp_norm_square_full["B2g"])
                d = -numpy.sqrt(symm_comp_norm_square_full["B3g"])
                e = numpy.sqrt(symm_comp_norm_square_full["Ag"])
                f = -numpy.sqrt(symm_comp_norm_square_full["B1g"])
            elif argmax_symm_comp_irrep_name in ["B2u", "B3u"]:
                c = numpy.sqrt(symm_comp_norm_square_full["B2u"])
                d = -numpy.sqrt(symm_comp_norm_square_full["B3u"])
                e = numpy.sqrt(symm_comp_norm_square_full["Au"])
                f = -numpy.sqrt(symm_comp_norm_square_full["B1u"])
            else:
                raise RuntimeError
        elif argmax_loc == 3:
            if argmax_symm_comp_irrep_name in ["B2g", "B3g"]:
                c = numpy.sqrt(symm_comp_norm_square_full["B3g"])
                d = numpy.sqrt(symm_comp_norm_square_full["B2g"])
                e = numpy.sqrt(symm_comp_norm_square_full["B1g"])
                f = numpy.sqrt(symm_comp_norm_square_full["Ag"])
            elif argmax_symm_comp_irrep_name in ["B2u", "B3u"]:
                c = numpy.sqrt(symm_comp_norm_square_full["B3u"])
                d = numpy.sqrt(symm_comp_norm_square_full["B2u"])
                e = numpy.sqrt(symm_comp_norm_square_full["B1u"])
                f = numpy.sqrt(symm_comp_norm_square_full["Au"])
            else:
                raise RuntimeError
        else:
            raise RuntimeError

        print(res[:, i])
        print(res[:, i + 1])

        new_basis1 = res[:, i] * (c - d * (1.0j)) + res[:, i + 1] * (-e - f * (1.0j))
        new_basis2 = res[:, i] * (e - f * (1.0j)) + res[:, i + 1] * (c + d * (1.0j))
        res[:, i], res[:, i + 1] = new_basis1, new_basis2

        print(res[:, i])
        print(res[:, i + 1])

        # print(c, d, e, f)
        print(
            "transform with factor = %20.16f + %20.16f i, %20.16f + %20.16f i"
            % (c, d, e, f)
        )

        mo_Lalpha_coeff = ao2Lalpha @ res[:n2c, :]
        mo_Lbeta_coeff = ao2Lbeta @ res[:n2c, :]
        mo_Lalpha_coeff_real = mo_Lalpha_coeff.real
        mo_Lalpha_coeff_imag = mo_Lalpha_coeff.imag
        mo_Lbeta_coeff_real = mo_Lbeta_coeff.real
        mo_Lbeta_coeff_imag = mo_Lbeta_coeff.imag

    # final check #

    n2c = mol.nao_2c()
    mo_Lalpha_coeff = ao2Lalpha @ res[:n2c, :]
    mo_Lbeta_coeff = ao2Lbeta @ res[:n2c, :]
    mo_Lalpha_coeff_real = mo_Lalpha_coeff.real
    mo_Lalpha_coeff_imag = mo_Lalpha_coeff.imag
    mo_Lbeta_coeff_real = mo_Lbeta_coeff.real
    mo_Lbeta_coeff_imag = mo_Lbeta_coeff.imag

    for i in range(0, n2c, 2):  # loop over unbarred only #

        print("---------------------------------------")
        print("i = ", i)

        # print(res[:, i])
        # print(res[:, i+1])

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

        ### judge whether it is purified ###

        # argmax_loc = numpy.argmax(numpy.array(norm_lst))

        print(norm_lst)

        for argmax_loc in range(4):
            symm_comp_norm_square, _ = label_orb_symm(
                mol,
                mol.irrep_name,
                mol.symm_orb,
                mo_coeff_lst[argmax_loc].reshape(-1, 1),
                ovlp,
                # tol=1e-7/norm_lst[argmax_loc],
            )
            argmax_symm_comp = numpy.argmax(numpy.array(symm_comp_norm_square))
            argmax_symm_comp_value = symm_comp_norm_square[argmax_symm_comp]
            argmax_symm_comp_irrep_name = mol.irrep_name[argmax_symm_comp]
            print(symm_comp_norm_square)
            print("argmax_loc ", argmax_loc)
            print(
                "symm comp %s %20.16f"
                % (argmax_symm_comp_irrep_name, argmax_symm_comp_value[0])
            )

        norm_alpha_real = (
            mo_Lalpha_coeff_real[:, i + 1].T @ ovlp @ mo_Lalpha_coeff_real[:, i + 1]
        )
        norm_alpha_imag = (
            mo_Lalpha_coeff_imag[:, i + 1].T @ ovlp @ mo_Lalpha_coeff_imag[:, i + 1]
        )
        norm_beta_real = (
            mo_Lbeta_coeff_real[:, i + 1].T @ ovlp @ mo_Lbeta_coeff_real[:, i + 1]
        )
        norm_beta_imag = (
            mo_Lbeta_coeff_imag[:, i + 1].T @ ovlp @ mo_Lbeta_coeff_imag[:, i + 1]
        )

        norm_lst = [norm_alpha_real, norm_alpha_imag, norm_beta_real, norm_beta_imag]
        mo_coeff_lst = numpy.array(
            [
                mo_Lalpha_coeff_real[:, i + 1],
                mo_Lalpha_coeff_imag[:, i + 1],
                mo_Lbeta_coeff_real[:, i + 1],
                mo_Lbeta_coeff_imag[:, i + 1],
            ]
        )

        print(norm_lst)

        ### judge whether it is purified ###

        # argmax_loc = numpy.argmax(numpy.array(norm_lst))

        for argmax_loc in range(4):

            symm_comp_norm_square, _ = label_orb_symm(
                mol,
                mol.irrep_name,
                mol.symm_orb,
                mo_coeff_lst[argmax_loc].reshape(-1, 1),
                ovlp,
                # tol=1e-7/norm_lst[argmax_loc],
            )

            argmax_symm_comp = numpy.argmax(numpy.array(symm_comp_norm_square))
            argmax_symm_comp_value = symm_comp_norm_square[argmax_symm_comp]
            argmax_symm_comp_irrep_name = mol.irrep_name[argmax_symm_comp]

            print(symm_comp_norm_square)

            print("argmax_loc ", argmax_loc)
            print(
                "symm comp %s %20.16f"
                % (argmax_symm_comp_irrep_name, argmax_symm_comp_value[0])
            )

    return res


def purify_pes_d2h_v2(mol, mf, pure_tol=1e-10, non_vanishing_tol=1e-8):

    IRREP_NAME_D2H_FULL = ["Ag", "B1g", "B2g", "B3g", "Au", "B1u", "B2u", "B3u"]

    ALLOWED_SYM = [["Ag", "Au"], ["B1g", "B1u"], ["B2g", "B2u"], ["B3g", "B3u"]]

    ALLOWED_SYM_UNBARRED = [
        ["B2g", "B2u"],
        ["B3g", "B3u"],
        ["Ag", "Au"],
        ["B1g", "B1u"],
    ]

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

    for i in range(0, n2c, 2):  # loop over unbarred only #

        print("---------------------------------------")
        print("i = ", i)

        #####################################
        ## collect the basic info ##
        #####################################

        while True:

            norm_alpha_real = (
                mo_Lalpha_coeff_real[:, i].T @ ovlp @ mo_Lalpha_coeff_real[:, i]
            )
            norm_alpha_imag = (
                mo_Lalpha_coeff_imag[:, i].T @ ovlp @ mo_Lalpha_coeff_imag[:, i]
            )
            norm_beta_real = (
                mo_Lbeta_coeff_real[:, i].T @ ovlp @ mo_Lbeta_coeff_real[:, i]
            )
            norm_beta_imag = (
                mo_Lbeta_coeff_imag[:, i].T @ ovlp @ mo_Lbeta_coeff_imag[:, i]
            )

            norm_lst = [
                norm_alpha_real,
                norm_alpha_imag,
                norm_beta_real,
                norm_beta_imag,
            ]
            mo_coeff_lst = numpy.array(
                [
                    mo_Lalpha_coeff_real[:, i],
                    mo_Lalpha_coeff_imag[:, i],
                    mo_Lbeta_coeff_real[:, i],
                    mo_Lbeta_coeff_imag[:, i],
                ]
            )

            norm_alpha_real_bar = (
                mo_Lalpha_coeff_real[:, i + 1].T @ ovlp @ mo_Lalpha_coeff_real[:, i + 1]
            )
            norm_alpha_imag_bar = (
                mo_Lalpha_coeff_imag[:, i + 1].T @ ovlp @ mo_Lalpha_coeff_imag[:, i + 1]
            )
            norm_beta_real_bar = (
                mo_Lbeta_coeff_real[:, i + 1].T @ ovlp @ mo_Lbeta_coeff_real[:, i + 1]
            )
            norm_beta_imag_bar = (
                mo_Lbeta_coeff_imag[:, i + 1].T @ ovlp @ mo_Lbeta_coeff_imag[:, i + 1]
            )

            norm_lst_bar = [
                norm_alpha_real_bar,
                norm_alpha_imag_bar,
                norm_beta_real_bar,
                norm_beta_imag_bar,
            ]
            mo_coeff_lst_bar = numpy.array(
                [
                    mo_Lalpha_coeff_real[:, i + 1],
                    mo_Lalpha_coeff_imag[:, i + 1],
                    mo_Lbeta_coeff_real[:, i + 1],
                    mo_Lbeta_coeff_imag[:, i + 1],
                ]
            )

            is_all_pure = True

            is_unbarred_pure = True
            is_barred_pure = True
            is_break_out = False

            for icomp in range(4):

                print("icomp = ", icomp)

                norm_comp = norm_lst[icomp]
                norm_comp_bar = norm_lst_bar[icomp]

                print("norm = ", norm_comp, norm_comp_bar)

                unbarred_non_vanishing = norm_comp > non_vanishing_tol
                bar_non_vanishing = norm_comp_bar > non_vanishing_tol

                if not unbarred_non_vanishing and not bar_non_vanishing:
                    continue

                if unbarred_non_vanishing and bar_non_vanishing:
                    symm_comp_norm_square, symorb = label_orb_symm(
                        mol,
                        mol.irrep_name,
                        mol.symm_orb,
                        mo_coeff_lst[icomp].reshape(-1, 1),
                        ovlp,
                    )

                    symm_comp_norm_square_bar, symorb_bar = label_orb_symm(
                        mol,
                        mol.irrep_name,
                        mol.symm_orb,
                        mo_coeff_lst_bar[icomp].reshape(-1, 1),
                        ovlp,
                    )

                    argmax_symm_comp = numpy.argmax(numpy.array(symm_comp_norm_square))
                    argmax_symm_comp_value = symm_comp_norm_square[argmax_symm_comp]
                    argmax_symm_comp_irrep_name = mol.irrep_name[argmax_symm_comp]

                    if argmax_symm_comp_value < (1 - pure_tol):
                        is_unbarred_pure = False

                    argmax_symm_comp_bar = numpy.argmax(
                        numpy.array(symm_comp_norm_square_bar)
                    )
                    argmax_symm_comp_value_bar = symm_comp_norm_square_bar[
                        argmax_symm_comp_bar
                    ]
                    argmax_symm_comp_irrep_name_bar = mol.irrep_name[
                        argmax_symm_comp_bar
                    ]

                    if argmax_symm_comp_value_bar < (1 - pure_tol):
                        is_barred_pure = False

                    if is_unbarred_pure and is_barred_pure:
                        print(
                            "both pure with ",
                            argmax_symm_comp_irrep_name,
                            argmax_symm_comp_irrep_name_bar,
                        )
                        
                        print(symorb[:, argmax_symm_comp])
                        print(symorb_bar[:, argmax_symm_comp_bar])
                    
                    # TODO: rotate 一下 #

                # never check others #

                is_all_pure = True
                is_unbarred_pure = True
                is_barred_pure = True
                break

                # if unbarred_not_vanishing #

                if unbarred_non_vanishing:
                    symm_comp_norm_square, _ = label_orb_symm(
                        mol,
                        mol.irrep_name,
                        mol.symm_orb,
                        mo_coeff_lst[icomp].reshape(-1, 1),
                        ovlp,
                    )

                    print("unbar ", symm_comp_norm_square)

                    argmax_symm_comp = numpy.argmax(numpy.array(symm_comp_norm_square))
                    argmax_symm_comp_value = symm_comp_norm_square[argmax_symm_comp]
                    argmax_symm_comp_irrep_name = mol.irrep_name[argmax_symm_comp]

                    if argmax_symm_comp_value < (1 - pure_tol):
                        is_unbarred_pure = False
                        is_barred_pure = False
                        is_all_pure = False

                        print("unbarred comp %d does not admit one irrep" % (icomp))

                        break

                    # 只有 1 个对称性 #

                    if argmax_symm_comp_irrep_name not in ALLOWED_SYM[icomp]:
                        # is_all_pure = False
                        is_break_out = True
                        begin_idx = (icomp // 2) * 2
                        is_phase_i = False
                        for idx in range(begin_idx, begin_idx + 2):
                            if (
                                argmax_symm_comp_irrep_name in ALLOWED_SYM[idx]
                            ):  # 差了 i 相位
                                res[:, i] *= 1.0j
                                res[:, i + 1] *= -1.0j
                                is_phase_i = True

                                print("unbarred comp %d change phase i" % (icomp))

                                break
                        # alpha, beta 反了 #
                        if not is_phase_i:
                            res[:, i], res[:, i + 1] = res[:, i + 1], res[:, i]
                            res[:, i + 1] *= -1.0
                            print("unbarred comp %d swap basis" % (icomp))

                if is_break_out:
                    break

                # if barred_not_vanishing #

                if bar_non_vanishing:
                    symm_comp_norm_square, _ = label_orb_symm(
                        mol,
                        mol.irrep_name,
                        mol.symm_orb,
                        mo_coeff_lst_bar[icomp].reshape(-1, 1),
                        ovlp,
                    )

                    print("bar ", symm_comp_norm_square)

                    argmax_symm_comp = numpy.argmax(numpy.array(symm_comp_norm_square))
                    argmax_symm_comp_value = symm_comp_norm_square[argmax_symm_comp]
                    argmax_symm_comp_irrep_name = mol.irrep_name[argmax_symm_comp]

                    if argmax_symm_comp_value < (1 - pure_tol):
                        is_unbarred_pure = False
                        is_barred_pure = False
                        is_all_pure = False

                        print("barred comp %d does not admit one irrep" % (icomp))

                        print(
                            unbarred_non_vanishing,
                            bar_non_vanishing,
                            is_unbarred_pure,
                            is_barred_pure,
                            is_all_pure,
                        )

                        break

                    # 只有 1 个对称性 #

                    if argmax_symm_comp_irrep_name not in ALLOWED_SYM_UNBARRED[icomp]:
                        # is_all_pure = False
                        is_break_out = True
                        begin_idx = (icomp // 2) * 2
                        is_phase_i = False
                        for idx in range(begin_idx, begin_idx + 2):
                            if (
                                argmax_symm_comp_irrep_name in ALLOWED_SYM_UNBARRED[idx]
                            ):  # 差了 i 相位
                                res[:, i + 1] *= 1.0j
                                res[:, i] *= -1.0j
                                is_phase_i = True
                                print("barred comp %d change phase i" % (icomp))
                                break
                        # alpha, beta 反了 #
                        print(
                            "argmax_symm_comp_irrep_name", argmax_symm_comp_irrep_name
                        )
                        if not is_phase_i:
                            res[:, i], res[:, i + 1] = res[:, i + 1], res[:, i]
                            res[:, i + 1] *= -1.0
                            print("barred comp %d swap basis" % (icomp))

                print(
                    unbarred_non_vanishing,
                    bar_non_vanishing,
                    is_unbarred_pure,
                    is_barred_pure,
                    is_all_pure,
                )

                if is_break_out:
                    break

                if not is_unbarred_pure and not is_barred_pure:
                    print("both not pure")
                    break

                if not unbarred_non_vanishing and not is_barred_pure:
                    print("bar phase rotation we have not deal with")
                    break

                if not bar_non_vanishing and not is_unbarred_pure:
                    print("unbar phase rotation we have not deal with")
                    break

                print(
                    unbarred_non_vanishing,
                    bar_non_vanishing,
                    is_unbarred_pure,
                    is_barred_pure,
                    is_all_pure,
                )

                # raise RuntimeError

            # break

            print("before check whether to break out while")

            if is_all_pure and not is_break_out:
                assert is_unbarred_pure and is_barred_pure
                break
            else:
                # update info #
                if not is_unbarred_pure and not is_barred_pure:
                    print("both unbarred and barred not pure")
                    break

                mo_Lalpha_coeff = ao2Lalpha @ res[:n2c, :]
                mo_Lbeta_coeff = ao2Lbeta @ res[:n2c, :]
                mo_Lalpha_coeff_real = mo_Lalpha_coeff.real
                mo_Lalpha_coeff_imag = mo_Lalpha_coeff.imag
                mo_Lbeta_coeff_real = mo_Lbeta_coeff.real
                mo_Lbeta_coeff_imag = mo_Lbeta_coeff.imag


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
        basis="3-21g",
        # basis="cc-pvdz-dk",
        verbose=5,
        charge=0,
        spin=0,
        symmetry="d2h",
        # unit='bohr'
    )

    # print(get_ao2L_spinor(mol))

    ao2Lalpha = get_ao2L_spinor(mol)
    ao2Lbeta = get_ao2L_spinor(mol, False)

    # do scf #

    mol.build()
    mf = scf.dhf.RDHF(mol)
    mf.conv_tol = 1e-12
    mf.kernel()
    purify_pes_d2h(mol, mf)

    purify_pes_d2h_v2(mol, mf)

    # print(mol.symm_orb)
    # print(mol.irrep_name)

    # exit(1)

    # learn symm module line by line #

#     sym_org = mol._symm_orig
#     sym_axes = mol._symm_axes
#     group = mol.groupname
#
#     print("sym_org", sym_org)
#     print("sym_axes", sym_axes)
#     print("group", group)

#     from pyscf.symm import geom, param, basis
#
#     prop_atoms = mol.format_atom(mol._atom, sym_org, sym_axes, "Bohr")
#     eql_atom_ids = geom.symm_identical_atoms(group, prop_atoms)
#
#     print(prop_atoms)  # 原子 + 坐标
#     print(eql_atom_ids)
#
#     ops = numpy.asarray([param.D2H_OPS[op] for op in param.OPERATOR_TABLE[group]])
#     chartab = numpy.array([x[1:] for x in param.CHARACTER_TABLE[group]])
#     nirrep = chartab.__len__()
#     aoslice = mol.aoslice_by_atom()
#     nao = mol.nao_nr()
#     atom_coords = numpy.array([a[1] for a in prop_atoms])
#
#     print(ops)  # 以 (x,y,z) 为 basis 的操作矩阵
#     print(chartab)
#     print(aoslice)
#     print(atom_coords)
#
#     print("------------------------------------------------")
#
#     sodic = [[] for i in range(8)]
#     for atom_ids in eql_atom_ids:
#         print("atom_ids = ", atom_ids)
#         r0 = atom_coords[atom_ids[0]]
#         print("r0 = ", r0)
#         op_coords = numpy.einsum(
#             "x,nxy->ny", r0, ops
#         )  # ops  作用到第一个原子得到的其他坐标
#         print(op_coords)
#         # Using ops to generate other atoms from atom_ids[0]
#         coords0 = atom_coords[atom_ids]
#         print("coords0 = ", coords0)
#         print(op_coords.shape, coords0.shape)
#         print((op_coords.reshape(-1, 1, 3) - coords0).shape)
#         dc = abs(op_coords.reshape(-1, 1, 3) - coords0).max(axis=2)
#         print("dc = ", dc)  # 记录操作关联了哪些原子
#         op_relate_idx = numpy.argwhere(dc < geom.TOLERANCE)[:, 1]
#         print("op_relate_idx = ", op_relate_idx)
#         ao_loc = numpy.array([aoslice[atom_ids[i], 2] for i in op_relate_idx])
#         print("ao_loc = ", ao_loc)  # 对称操作关联的atom 对应的basis开始的指标
#
#         b0, b1 = aoslice[atom_ids[0], :2]
#         ip = 0
#         for ib in range(b0, b1):
#             l = mol.bas_angular(ib)
#             print("ib = ", ib, "l = ", l)
#             if mol.cart:  # for spinor 报错 !
#                 degen = (l + 1) * (l + 2) // 2
#                 cbase = numpy.zeros((degen, nirrep, nao))
#                 for op_id, op in enumerate(ops):
#                     n = 0
#                     for x in range(l, -1, -1):
#                         for y in range(l - x, -1, -1):
#                             z = l - x - y
#                             idx = ao_loc[op_id] + n
#                             sign = op[0, 0] ** x * op[1, 1] ** y * op[2, 2] ** z
#                             cbase[n, :, idx] += sign * chartab[:, op_id]
#                             n += 1
#             else:
#                 degen = l * 2 + 1
#                 cbase = numpy.zeros((degen, nirrep, nao))
#                 for op_id, op in enumerate(param.OPERATOR_TABLE[group]):
#                     print("op_id = ", op_id, "op = ", op)
#                     for n, m in enumerate(range(-l, l + 1)):
#                         idx = ao_loc[op_id] + n
#                         print("n = ", n, "m = ", m, "idx = ", idx)
#                         if basis.tot_parity_odd(op, l, m):
#                             cbase[n, :, idx] -= chartab[:, op_id]
#                         else:
#                             cbase[n, :, idx] += chartab[:, op_id]
#             norms = numpy.sqrt(numpy.einsum("mij,mij->mi", cbase, cbase))
#
#             for i in range(mol.bas_nctr(ib)):
#                 for n, ir in numpy.argwhere(norms > 1e-12):
#                     c = numpy.zeros(nao)
#                     c[ip:] = cbase[n, ir, : nao - ip] / norms[n, ir]
#                     sodic[ir].append(c)
#                 ip += degen
#
#         print("------------------------------------------------")
#

# Step I 测试纯化程序正确定 #

# Step II 产生随机扰动检查程序是否正确 #
