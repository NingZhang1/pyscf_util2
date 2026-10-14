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


########################################
# 检查 spinor 是否 adapted #
########################################


def check_spinor_adapted(mol, pes_mocoeff, pure_tol=1e-10, non_vanishing_tol=1e-8):

    print(" ---------------------- check_spinor_adapted ---------------------- ")

    ALLOWED_SYM = [["Ag", "Au"], ["B1g", "B1u"], ["B2g", "B2u"], ["B3g", "B3u"]]
    ALLOWED_SYM_BARRED = [
        ["B2g", "B2u"],
        ["B3g", "B3u"],
        ["Ag", "Au"],
        ["B1g", "B1u"],
    ]

    # generate ao2L component #

    ao2Lalpha = get_ao2L_spinor(mol)
    ao2Lbeta = get_ao2L_spinor(mol, False)

    # generate mo_Lalpha, mo_Lbeta #

    n2c = pes_mocoeff.shape[1]
    mo_Lalpha_coeff = ao2Lalpha @ pes_mocoeff[:n2c, :]
    mo_Lbeta_coeff = ao2Lbeta @ pes_mocoeff[:n2c, :]
    mo_Lalpha_coeff_real = mo_Lalpha_coeff.real
    mo_Lalpha_coeff_imag = mo_Lalpha_coeff.imag
    mo_Lbeta_coeff_real = mo_Lbeta_coeff.real
    mo_Lbeta_coeff_imag = mo_Lbeta_coeff.imag

    # res = pes_mocoeff.copy()

    ovlp = mol.intor("int1e_ovlp")
    result = []
    result_dominant_comp = []
    for _ in range(n2c // 2):
        result_dominant_comp.append(1.0)

    print("mol.irrep_name = ", mol.irrep_name)

    for i in range(0, n2c, 2):

        print("do i = ", i)

        norm_alpha_real = (
            mo_Lalpha_coeff_real[:, i].T @ ovlp @ mo_Lalpha_coeff_real[:, i]
        )
        norm_alpha_imag = (
            mo_Lalpha_coeff_imag[:, i].T @ ovlp @ mo_Lalpha_coeff_imag[:, i]
        )
        norm_beta_real = mo_Lbeta_coeff_real[:, i].T @ ovlp @ mo_Lbeta_coeff_real[:, i]
        norm_beta_imag = mo_Lbeta_coeff_imag[:, i].T @ ovlp @ mo_Lbeta_coeff_imag[:, i]

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

        is_all_adapted = True

        # print("norm_lst     = ", norm_lst)
        # print("norm_lst_bar = ", norm_lst_bar)

        argmax_loc = numpy.argmax(numpy.array(norm_lst))
        argmax_loc_bar = numpy.argmax(numpy.array(norm_lst_bar))

        unbar_not_pure_norm = None
        bar_not_pure_norm = None
        argmax_symm_comp_value_not_pure = None
        argmax_symm_comp_value_bar_not_pure = None

        for icomp in range(4):

            # print("icomp = ", icomp)

            norm_comp = norm_lst[icomp]
            norm_comp_bar = norm_lst_bar[icomp]

            unbarred_non_vanishing = norm_comp > non_vanishing_tol
            bar_non_vanishing = norm_comp_bar > non_vanishing_tol

            if not unbarred_non_vanishing and not bar_non_vanishing:
                continue

            if unbarred_non_vanishing:
                symm_comp_norm_square, _ = label_orb_symm(
                    mol,
                    mol.irrep_name,
                    mol.symm_orb,
                    mo_coeff_lst[icomp].reshape(-1, 1),
                    ovlp,
                )

                argmax_symm_comp = numpy.argmax(numpy.array(symm_comp_norm_square))
                argmax_symm_comp_value = symm_comp_norm_square[argmax_symm_comp]
                argmax_symm_comp_irrep_name = mol.irrep_name[argmax_symm_comp]

                # print("unbarred = ", symm_comp_norm_square)

                if argmax_symm_comp_value < (1 - pure_tol):
                    print(
                        "comp %d unbarred is not pure since no single irrep dominates with largest comp %20.16f"
                        % (icomp, argmax_symm_comp_value)
                    )
                    is_all_adapted = False

                    # result_dominant_comp[i // 2] = argmax_symm_comp_value

                    if unbar_not_pure_norm is None:
                        unbar_not_pure_norm = norm_lst[icomp]
                        result_dominant_comp[i // 2] = argmax_symm_comp_value
                        argmax_symm_comp_value_not_pure = argmax_symm_comp_value
                    else:
                        if norm_lst[icomp] > unbar_not_pure_norm:
                            unbar_not_pure_norm = norm_lst[icomp]
                            result_dominant_comp[i // 2] = argmax_symm_comp_value
                            argmax_symm_comp_value_not_pure = argmax_symm_comp_value

                if is_all_adapted:
                    if argmax_symm_comp_irrep_name not in ALLOWED_SYM[icomp]:
                        print(
                            "comp %d unbarred is not pure since the irrep %s is not in %s"
                            % (icomp, argmax_symm_comp_irrep_name, ALLOWED_SYM[icomp])
                        )
                        is_all_adapted = False
                        result_dominant_comp[i // 2] = 0.0

            # if not is_all_adapted:
            #     break

            if bar_non_vanishing:
                symm_comp_norm_square_bar, _ = label_orb_symm(
                    mol,
                    mol.irrep_name,
                    mol.symm_orb,
                    mo_coeff_lst_bar[icomp].reshape(-1, 1),
                    ovlp,
                )

                # print("barred = ", symm_comp_norm_square_bar)

                argmax_symm_comp = numpy.argmax(numpy.array(symm_comp_norm_square_bar))
                argmax_symm_comp_value = symm_comp_norm_square_bar[argmax_symm_comp]
                argmax_symm_comp_irrep_name = mol.irrep_name[argmax_symm_comp]

                if argmax_symm_comp_value < (1 - pure_tol):
                    print(
                        "comp %d barred is not pure since no single irrep dominates with largest comp %20.16f"
                        % (icomp, argmax_symm_comp_value)
                    )
                    is_all_adapted = False
                    # result_dominant_comp[i // 2] = argmax_symm_comp_value

                    if bar_not_pure_norm is None:
                        bar_not_pure_norm = norm_lst_bar[icomp]
                        result_dominant_comp[i // 2] = argmax_symm_comp_value
                        argmax_symm_comp_value_bar_not_pure = argmax_symm_comp_value
                    else:
                        if norm_lst_bar[icomp] > bar_not_pure_norm:
                            bar_not_pure_norm = norm_lst_bar[icomp]
                            result_dominant_comp[i // 2] = argmax_symm_comp_value
                            argmax_symm_comp_value_bar_not_pure = argmax_symm_comp_value

                if is_all_adapted:
                    if argmax_symm_comp_irrep_name not in ALLOWED_SYM_BARRED[icomp]:
                        print(
                            "comp %d barred is not pure since the irrep %s is not in %s"
                            % (
                                icomp,
                                argmax_symm_comp_irrep_name,
                                ALLOWED_SYM_BARRED[icomp],
                            )
                        )
                        # print(norm_lst)
                        # print(norm_lst_bar)
                        # print(symm_comp_norm_square)
                        # print(symm_comp_norm_square_bar)
                        is_all_adapted = False
                        result_dominant_comp[i // 2] = 0.0

            # if not is_all_adapted:
            #     break

            if (
                argmax_symm_comp_value_bar_not_pure is not None
                or argmax_symm_comp_value_not_pure is not None
            ):
                assert not is_all_adapted
                if (
                    argmax_symm_comp_value_bar_not_pure is not None
                    and argmax_symm_comp_value_not_pure is not None
                ):
                    result_dominant_comp[i // 2] = max(
                        argmax_symm_comp_value_bar_not_pure,
                        argmax_symm_comp_value_not_pure,
                    )
                else:
                    if argmax_symm_comp_value_bar_not_pure is not None:
                        result_dominant_comp[i // 2] = (
                            argmax_symm_comp_value_bar_not_pure
                        )
                    else:
                        result_dominant_comp[i // 2] = argmax_symm_comp_value_not_pure

        result.append(is_all_adapted)

    print(" ---------------------- check_spinor_adapted end ---------------------- ")

    return result, result_dominant_comp


########################################
# 处理 symmetry pure 情形 #
########################################


def purify_spinor_with_dominant_irrep(
    mol, pes_mocoeff, adapted_spinor, pure_tol=1e-10, non_vanishing_tol=1e-8
):

    # nothing to do #

    if all(adapted_spinor):
        return pes_mocoeff

    print(
        " ---------------------- purify_spinor_with_dominant_irrep ---------------------- "
    )

    ALLOWED_SYM = [["Ag", "Au"], ["B1g", "B1u"], ["B2g", "B2u"], ["B3g", "B3u"]]
    ALLOWED_SYM_BARRED = [
        ["B2g", "B2u"],
        ["B3g", "B3u"],
        ["Ag", "Au"],
        ["B1g", "B1u"],
    ]

    # generate ao2L component #

    ao2Lalpha = get_ao2L_spinor(mol)
    ao2Lbeta = get_ao2L_spinor(mol, False)

    # generate mo_Lalpha, mo_Lbeta #

    n2c = pes_mocoeff.shape[1]
    mo_Lalpha_coeff = ao2Lalpha @ pes_mocoeff[:n2c, :]
    mo_Lbeta_coeff = ao2Lbeta @ pes_mocoeff[:n2c, :]
    mo_Lalpha_coeff_real = mo_Lalpha_coeff.real
    mo_Lalpha_coeff_imag = mo_Lalpha_coeff.imag
    mo_Lbeta_coeff_real = mo_Lbeta_coeff.real
    mo_Lbeta_coeff_imag = mo_Lbeta_coeff.imag

    res = pes_mocoeff.copy()

    ovlp = mol.intor("int1e_ovlp")

    for i in range(0, n2c, 2):

        if adapted_spinor[i // 2]:
            continue

        print("do i = ", i)

        can_be_fixed = True  # by x phase i or swap alpha beta

        norm_alpha_real = (
            mo_Lalpha_coeff_real[:, i].T @ ovlp @ mo_Lalpha_coeff_real[:, i]
        )
        norm_alpha_imag = (
            mo_Lalpha_coeff_imag[:, i].T @ ovlp @ mo_Lalpha_coeff_imag[:, i]
        )
        norm_beta_real = mo_Lbeta_coeff_real[:, i].T @ ovlp @ mo_Lbeta_coeff_real[:, i]
        norm_beta_imag = mo_Lbeta_coeff_imag[:, i].T @ ovlp @ mo_Lbeta_coeff_imag[:, i]

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

        for icomp in range(4):

            # print("icomp = ", icomp)

            norm_comp = norm_lst[icomp]
            norm_comp_bar = norm_lst_bar[icomp]

            # print("norm = ", norm_comp, norm_comp_bar)

            unbarred_non_vanishing = norm_comp > non_vanishing_tol
            bar_non_vanishing = norm_comp_bar > non_vanishing_tol

            if not unbarred_non_vanishing and not bar_non_vanishing:
                continue

            if unbarred_non_vanishing:
                symm_comp_norm_square, _ = label_orb_symm(
                    mol,
                    mol.irrep_name,
                    mol.symm_orb,
                    mo_coeff_lst[icomp].reshape(-1, 1),
                    ovlp,
                )

                argmax_symm_comp = numpy.argmax(numpy.array(symm_comp_norm_square))
                argmax_symm_comp_value = symm_comp_norm_square[argmax_symm_comp]
                argmax_symm_comp_irrep_name = mol.irrep_name[argmax_symm_comp]

                if argmax_symm_comp_value < (1 - pure_tol):
                    print(
                        "comp %d unbarred can not be fixed since no single irrep dominates"
                        % icomp
                    )
                    can_be_fixed = False

                # if can_be_fixed:
                #     if argmax_symm_comp_irrep_name not in ALLOWED_SYM[icomp]:
                #         begin_idx = (icomp // 2) * 2
                #         found = False
                #         for idx in range(begin_idx, begin_idx + 2):
                #             if argmax_symm_comp_irrep_name not in ALLOWED_SYM[idx]:
                #                 found = True
                #                 break
                #         if not found:
                #             can_be_fixed = False

            if not can_be_fixed:
                break

            if bar_non_vanishing:
                symm_comp_norm_square, _ = label_orb_symm(
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
                    print(
                        "comp %d barred can not be fixed since no single irrep dominates"
                        % icomp
                    )
                    can_be_fixed = False

                # if can_be_fixed:
                #     if argmax_symm_comp_irrep_name not in ALLOWED_SYM_BARRED[icomp]:
                #         begin_idx = (icomp // 2) * 2
                #         found = False
                #         for idx in range(begin_idx, begin_idx + 2):
                #             if (
                #                 argmax_symm_comp_irrep_name
                #                 not in ALLOWED_SYM_BARRED[idx]
                #             ):
                #                 found = True
                #                 break
                #         if not found:
                #             can_be_fixed = False

            if not can_be_fixed:
                break

        if can_be_fixed:

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
            # argmax_symm_comp_value = symm_comp_norm_square[argmax_symm_comp]
            argmax_symm_comp_irrep_name = mol.irrep_name[argmax_symm_comp]

            norm_comp = norm_lst[argmax_loc]
            norm_comp_bar = norm_lst_bar[argmax_loc]

            if norm_comp_bar > non_vanishing_tol:
                symm_comp_norm_square_bar, _ = label_orb_symm(
                    mol,
                    mol.irrep_name,
                    mol.symm_orb,
                    mo_coeff_lst_bar[argmax_loc].reshape(-1, 1),
                    ovlp,
                    # tol=1e-7/norm_lst[argmax_loc],
                )

                argmax_symm_comp_bar = numpy.argmax(
                    numpy.array(symm_comp_norm_square_bar)
                )
                # argmax_symm_comp_value_bar = symm_comp_norm_square_bar[argmax_symm_comp_bar]
                argmax_symm_comp_irrep_name_bar = mol.irrep_name[argmax_symm_comp_bar]

            # first judge whether should be multiply by phase i #

            # assert argmax_symm_comp_irrep_name not in ALLOWED_SYM[argmax_loc]

            TO_OTHER = [1, 0, 3, 2]
            argmax_loc_other = TO_OTHER[argmax_loc]
            TO_SWAP = [2, 3, 0, 1]
            TO_OTHER_SWAP = [3, 2, 1, 0]

            if (
                argmax_symm_comp_irrep_name in ALLOWED_SYM[argmax_loc_other]
            ):  # multiby i
                can_be_fixed_by_multiply_i = True
                if norm_comp_bar > non_vanishing_tol:
                    if (
                        argmax_symm_comp_irrep_name_bar
                        not in ALLOWED_SYM_BARRED[argmax_loc_other]
                    ):
                        can_be_fixed_by_multiply_i = False
                        can_be_fixed = False

                if can_be_fixed_by_multiply_i:
                    res[:, i] *= 1.0j
                    res[:, i + 1] *= -1.0j
                    print("fixed by multiply by i")
                    continue
                else:
                    print("cannot be fixed by multiply by i")

            if not can_be_fixed:
                continue

            if (
                argmax_symm_comp_irrep_name in ALLOWED_SYM[TO_SWAP[argmax_loc]]
            ):  # swap alpha beta
                can_be_fixed_by_swap = True
                if norm_comp_bar > non_vanishing_tol:
                    if (
                        argmax_symm_comp_irrep_name_bar
                        not in ALLOWED_SYM_BARRED[TO_SWAP[argmax_loc]]
                    ):
                        can_be_fixed_by_swap = False
                        can_be_fixed = False

                if can_be_fixed_by_swap:
                    res[:, i], res[:, i + 1] = res[:, i + 1], res[:, i]
                    res[:, i + 1] *= -1.0
                    print("swap alpha beta")
                    continue
                else:
                    print("cannot be fixed by swap alpha beta")

            if not can_be_fixed:
                continue

            if (
                argmax_symm_comp_irrep_name in ALLOWED_SYM[TO_OTHER_SWAP[argmax_loc]]
            ):  # swap alpha beta and multiply by i
                can_be_fixed_by_swap_i = True
                if norm_comp_bar > non_vanishing_tol:
                    if (
                        argmax_symm_comp_irrep_name_bar
                        not in ALLOWED_SYM_BARRED[TO_OTHER_SWAP[argmax_loc]]
                    ):
                        can_be_fixed_by_swap_i = False
                        can_be_fixed = False

                if can_be_fixed_by_swap_i:
                    res[:, i], res[:, i + 1] = res[:, i + 1], res[:, i]
                    res[:, i + 1] *= -1.0
                    res[:, i] *= 1.0j
                    res[:, i + 1] *= -1.0j
                    print("swap alpha beta and multiply by i")
                    continue
                else:
                    print("cannot be fixed by swap alpha beta and multiply by i")

    print(
        " ---------------------- purify_spinor_with_dominant_irrep end ---------------------- "
    )

    return res


########################################
# 处理 需要混合的 情形 #
# symmetry pure 情形 其实是这个的特例 #
########################################


def purify_spinor_with_mixed_irrep(
    mol,
    pes_mocoeff,
    adapted_spinor,
    pure_tol=1e-10,
    non_vanishing_tol=1e-8,
    always_use_max_norm_comp=False,
):

    # nothing to do #

    if all(adapted_spinor):
        return pes_mocoeff

    print(
        " ---------------------- purify_spinor_with_mixed_irrep ---------------------- "
    )

    IRREP_NAME_G = ["Ag", "B1g", "B2g", "B3g"]
    IRREP_NAME_D2H_FULL = ["Ag", "B1g", "B2g", "B3g", "Au", "B1u", "B2u", "B3u"]

    # generate ao2L component #

    ao2Lalpha = get_ao2L_spinor(mol)
    ao2Lbeta = get_ao2L_spinor(mol, False)

    # generate mo_Lalpha, mo_Lbeta #

    n2c = pes_mocoeff.shape[1]
    mo_Lalpha_coeff = ao2Lalpha @ pes_mocoeff[:n2c, :]
    mo_Lbeta_coeff = ao2Lbeta @ pes_mocoeff[:n2c, :]
    mo_Lalpha_coeff_real = mo_Lalpha_coeff.real
    mo_Lalpha_coeff_imag = mo_Lalpha_coeff.imag
    mo_Lbeta_coeff_real = mo_Lbeta_coeff.real
    mo_Lbeta_coeff_imag = mo_Lbeta_coeff.imag

    res = pes_mocoeff.copy()
    ovlp = mol.intor("int1e_ovlp")

    for i in range(0, n2c, 2):

        if adapted_spinor[i // 2]:
            continue

        print("------------------------------------------------")
        print("do i = ", i)

        norm_alpha_real = (
            mo_Lalpha_coeff_real[:, i].T @ ovlp @ mo_Lalpha_coeff_real[:, i]
        )
        norm_alpha_imag = (
            mo_Lalpha_coeff_imag[:, i].T @ ovlp @ mo_Lalpha_coeff_imag[:, i]
        )
        norm_beta_real = mo_Lbeta_coeff_real[:, i].T @ ovlp @ mo_Lbeta_coeff_real[:, i]
        norm_beta_imag = mo_Lbeta_coeff_imag[:, i].T @ ovlp @ mo_Lbeta_coeff_imag[:, i]

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

        # print(norm_lst)
        # print(norm_lst_bar)

        norm_test = [norm_lst[0], norm_lst[1], norm_lst_bar[0], norm_lst_bar[1]]
        if all([x < non_vanishing_tol for x in norm_test]):
            raise ValueError("Impossible Input")
        # print(norm_test)

        IsG = None
        # symorb = None
        symorb = {}
        for x in IRREP_NAME_D2H_FULL:
            symorb[x] = None

        if norm_test[0] > non_vanishing_tol:

            symm_comp_norm_square_unbar_real, symmorb_unbar_real = label_orb_symm(
                mol,
                mol.irrep_name,
                mol.symm_orb,
                mo_coeff_lst[0].reshape(-1, 1),
                ovlp,
            )

            symm_comp_norm_square_unbar_real_full = {}
            symmorb_unbar_real_full = {}
            for x, name in enumerate(mol.irrep_name):
                symm_comp_norm_square_unbar_real_full[name] = (
                    symm_comp_norm_square_unbar_real[x]
                )
                symmorb_unbar_real_full[name] = symmorb_unbar_real[:, x]

                if (
                    symm_comp_norm_square_unbar_real[x] * norm_test[0]
                ) > non_vanishing_tol and symorb[name] is None:
                    print("set symorb from unbar real ", name)
                    symorb[name] = symmorb_unbar_real[:, x]
            for x in IRREP_NAME_D2H_FULL:
                if x not in symm_comp_norm_square_unbar_real_full:
                    symm_comp_norm_square_unbar_real_full[x] = 0.0
                    symmorb_unbar_real_full[x] = None

            sum_g = 0.0
            for x in IRREP_NAME_G:
                sum_g += symm_comp_norm_square_unbar_real_full[x]
            if sum_g > (1.0 - pure_tol):
                IsG = True
            else:
                IsG = False

            # symorb = symmorb_unbar_real_full

        else:
            symm_comp_norm_square_unbar_real_full = {}
            for x in IRREP_NAME_D2H_FULL:
                symm_comp_norm_square_unbar_real_full[x] = 0.0
            symmorb_unbar_real_full = None

        if norm_test[1] > non_vanishing_tol:

            symm_comp_norm_square_unbar_imag, symmorb_unbar_imag = label_orb_symm(
                mol,
                mol.irrep_name,
                mol.symm_orb,
                mo_coeff_lst[1].reshape(-1, 1),
                ovlp,
            )

            # print(symmorb_unbar_imag)

            symm_comp_norm_square_unbar_imag_full = {}
            symmorb_unbar_imag_full = {}
            for x, name in enumerate(mol.irrep_name):
                symm_comp_norm_square_unbar_imag_full[name] = (
                    symm_comp_norm_square_unbar_imag[x]
                )
                symmorb_unbar_imag_full[name] = symmorb_unbar_imag[:, x]
                if (
                    symm_comp_norm_square_unbar_imag[x] * norm_test[1]
                ) > non_vanishing_tol and symorb[name] is None:
                    print("set symorb from unbar imag ", name)
                    symorb[name] = symmorb_unbar_imag[:, x]
            for x in IRREP_NAME_D2H_FULL:
                if x not in symm_comp_norm_square_unbar_imag_full:
                    symm_comp_norm_square_unbar_imag_full[x] = 0.0
                    symmorb_unbar_imag_full[x] = None

            sum_g = 0.0
            for x in IRREP_NAME_G:
                sum_g += symm_comp_norm_square_unbar_imag_full[x]
            if IsG is not None:
                if IsG:
                    assert sum_g > (1.0 - pure_tol)
                else:
                    # print(sum_g)
                    assert sum_g < pure_tol
            else:
                if sum_g > (1.0 - pure_tol):
                    IsG = True
                else:
                    IsG = False
                # symorb = symmorb_unbar_imag_full

        else:

            symm_comp_norm_square_unbar_imag_full = {}
            for x in IRREP_NAME_D2H_FULL:
                symm_comp_norm_square_unbar_imag_full[x] = 0.0
            symmorb_unbar_imag_full = None

        if norm_test[2] > non_vanishing_tol:

            symm_comp_norm_square_bar_real, symmorb_bar_real = label_orb_symm(
                mol,
                mol.irrep_name,
                mol.symm_orb,
                mo_coeff_lst_bar[0].reshape(-1, 1),
                ovlp,
            )

            symm_comp_norm_square_bar_real_full = {}
            symmorb_bar_real_full = {}
            for x, name in enumerate(mol.irrep_name):
                symm_comp_norm_square_bar_real_full[name] = (
                    symm_comp_norm_square_bar_real[x]
                )
                symmorb_bar_real_full[name] = symmorb_bar_real[:, x]
                if (
                    symm_comp_norm_square_bar_real[x] * norm_test[2]
                ) > non_vanishing_tol and symorb[name] is None:
                    print("set symorb from bar real ", name)
                    symorb[name] = symmorb_bar_real[:, x]
            for x in IRREP_NAME_D2H_FULL:
                if x not in symm_comp_norm_square_bar_real_full:
                    symm_comp_norm_square_bar_real_full[x] = 0.0
                    symmorb_bar_real_full[x] = None

            sum_g = 0.0
            for x in IRREP_NAME_G:
                sum_g += symm_comp_norm_square_bar_real_full[x]
            if IsG is not None:
                if IsG:
                    assert sum_g > (1.0 - pure_tol)
                else:
                    # print(sum_g)
                    assert sum_g < pure_tol
            else:
                if sum_g > (1.0 - pure_tol):
                    IsG = True
                else:
                    IsG = False
                # symorb = symmorb_bar_real_full

        else:

            symm_comp_norm_square_bar_real_full = {}
            for x in IRREP_NAME_D2H_FULL:
                symm_comp_norm_square_bar_real_full[x] = 0.0
            symmorb_bar_real_full = None

        if norm_test[3] > non_vanishing_tol:

            symm_comp_norm_square_bar_imag, symmorb_bar_imag = label_orb_symm(
                mol,
                mol.irrep_name,
                mol.symm_orb,
                mo_coeff_lst_bar[1].reshape(-1, 1),
                ovlp,
            )

            symm_comp_norm_square_bar_imag_full = {}
            symmorb_bar_imag_full = {}
            for x, name in enumerate(mol.irrep_name):
                symm_comp_norm_square_bar_imag_full[name] = (
                    symm_comp_norm_square_bar_imag[x]
                )

                symmorb_bar_imag_full[name] = symmorb_bar_imag[:, x]
                if (
                    symm_comp_norm_square_bar_imag[x] * norm_test[3]
                ) > non_vanishing_tol and symorb[name] is None:
                    print("set symorb from bar imag ", name)
                    symorb[name] = symmorb_bar_imag[:, x]
            for x in IRREP_NAME_D2H_FULL:
                if x not in symm_comp_norm_square_bar_imag_full:
                    symm_comp_norm_square_bar_imag_full[x] = 0.0
                    symmorb_bar_imag_full[x] = None

            sum_g = 0.0
            for x in IRREP_NAME_G:
                sum_g += symm_comp_norm_square_bar_imag_full[x]
            if IsG is not None:
                if IsG:
                    assert sum_g > (1.0 - pure_tol)
                else:
                    # print(sum_g)
                    assert sum_g < pure_tol
            else:
                if sum_g > (1.0 - pure_tol):
                    IsG = True
                else:
                    IsG = False
                # symorb = symmorb_bar_imag_full
        else:

            symm_comp_norm_square_bar_imag_full = {}
            for x in IRREP_NAME_D2H_FULL:
                symm_comp_norm_square_bar_imag_full[x] = 0.0
            symmorb_bar_imag_full = None

        assert IsG is not None

        for name in IRREP_NAME_D2H_FULL:
            if not isinstance(symm_comp_norm_square_unbar_real_full[name], float):
                symm_comp_norm_square_unbar_real_full[name] = (
                    symm_comp_norm_square_unbar_real_full[name][0]
                )
            if not isinstance(symm_comp_norm_square_unbar_imag_full[name], float):
                symm_comp_norm_square_unbar_imag_full[name] = (
                    symm_comp_norm_square_unbar_imag_full[name][0]
                )
            if not isinstance(symm_comp_norm_square_bar_real_full[name], float):
                symm_comp_norm_square_bar_real_full[name] = (
                    symm_comp_norm_square_bar_real_full[name][0]
                )
            if not isinstance(symm_comp_norm_square_bar_imag_full[name], float):
                symm_comp_norm_square_bar_imag_full[name] = (
                    symm_comp_norm_square_bar_imag_full[name][0]
                )

        #########################################
        # align sym orb #

        if symmorb_unbar_real_full is not None:

            for x in IRREP_NAME_D2H_FULL:
                if symorb[x] is not None and symmorb_unbar_real_full[x] is not None:
                    if symm_comp_norm_square_unbar_real_full[x] < non_vanishing_tol:
                        continue
                    print(
                        "unbar real ",
                        x,
                        norm_test[0],
                        symm_comp_norm_square_unbar_real_full[x],
                    )

                    sum_mat = symorb[x] - symmorb_unbar_real_full[x]
                    minus_mat = symorb[x] + symmorb_unbar_real_full[x]

                    max_norm_sum = numpy.max(numpy.abs(sum_mat))
                    max_norm_minus = numpy.max(numpy.abs(minus_mat))

                    if min(max_norm_sum, max_norm_minus) > non_vanishing_tol:
                        print(
                            "Warninig symorb not close with max_diff %10.4e"
                            % (min(max_norm_sum, max_norm_minus))
                        )

        if symmorb_unbar_imag_full is not None:
            for x in IRREP_NAME_D2H_FULL:
                if symorb[x] is not None and symmorb_unbar_imag_full[x] is not None:
                    if symm_comp_norm_square_unbar_imag_full[x] < non_vanishing_tol:
                        continue
                    print(
                        "unbar imag ",
                        x,
                        norm_test[1],
                        symm_comp_norm_square_unbar_imag_full[x],
                    )

                    sum_mat = symorb[x] - symmorb_unbar_imag_full[x]
                    minus_mat = symorb[x] + symmorb_unbar_imag_full[x]

                    max_norm_sum = numpy.max(numpy.abs(sum_mat))
                    max_norm_minus = numpy.max(numpy.abs(minus_mat))

                    if min(max_norm_sum, max_norm_minus) > non_vanishing_tol:
                        print(
                            "Warninig symorb not close with max_diff %10.4e"
                            % (min(max_norm_sum, max_norm_minus))
                        )

        if symmorb_bar_real_full is not None:
            for x in IRREP_NAME_D2H_FULL:
                if symorb[x] is not None and symmorb_bar_real_full[x] is not None:
                    if symm_comp_norm_square_bar_real_full[x] < non_vanishing_tol:
                        continue
                    print(
                        "bar real ",
                        x,
                        norm_test[2],
                        symm_comp_norm_square_bar_real_full[x],
                    )

                    sum_mat = symorb[x] - symmorb_bar_real_full[x]
                    minus_mat = symorb[x] + symmorb_bar_real_full[x]

                    max_norm_sum = numpy.max(numpy.abs(sum_mat))
                    max_norm_minus = numpy.max(numpy.abs(minus_mat))

                    if min(max_norm_sum, max_norm_minus) > non_vanishing_tol:
                        print(
                            "Warninig symorb not close with max_diff %10.4e"
                            % (min(max_norm_sum, max_norm_minus))
                        )

        if symmorb_bar_imag_full is not None:
            for x in IRREP_NAME_D2H_FULL:
                if symorb[x] is not None and symmorb_bar_imag_full[x] is not None:
                    if symm_comp_norm_square_bar_imag_full[x] < non_vanishing_tol:
                        continue
                    print(
                        "bar imag ",
                        x,
                        norm_test[3],
                        symm_comp_norm_square_bar_imag_full[x],
                    )

                    sum_mat = symorb[x] - symmorb_bar_imag_full[x]
                    minus_mat = symorb[x] + symmorb_bar_imag_full[x]

                    max_norm_sum = numpy.max(numpy.abs(sum_mat))
                    max_norm_minus = numpy.max(numpy.abs(minus_mat))

                    if min(max_norm_sum, max_norm_minus) > non_vanishing_tol:
                        print(
                            "Warninig symorb not close with max_diff %10.4e"
                            % (min(max_norm_sum, max_norm_minus))
                        )

        # exit(1)

        # print(symm_comp_norm_square_unbar_real_full)
        # print(symm_comp_norm_square_unbar_imag_full)
        # print(symm_comp_norm_square_bar_real_full)
        # print(symm_comp_norm_square_bar_imag_full)
        # print(norm_test)

        c = 0.0
        d = 0.0
        e = 0.0
        f = 0.0

        entropy = [0.0, 0.0, 0.0, 0.0]
        entropy = numpy.asarray(entropy)

        for x in IRREP_NAME_D2H_FULL:
            value = symm_comp_norm_square_unbar_real_full[x]
            if value > 1e-14:
                entropy[0] -= float(value * numpy.log(value))
            value = symm_comp_norm_square_unbar_imag_full[x]
            if value > 1e-14:
                entropy[1] -= float(value * numpy.log(value))
            value = symm_comp_norm_square_bar_real_full[x]
            if value > 1e-14:
                entropy[2] -= float(value * numpy.log(value))
            value = symm_comp_norm_square_bar_imag_full[x]
            if value > 1e-14:
                entropy[3] -= float(value * numpy.log(value))
        # entropy = numpy.asarray(entropy)
        # print("entropy = ", entropy)

        if IsG:

            norm_square_Ag = (
                norm_test[0] * symm_comp_norm_square_unbar_real_full["Ag"]
                + norm_test[1] * symm_comp_norm_square_unbar_imag_full["Ag"]
            )
            norm_square_Ag += (
                norm_test[2] * symm_comp_norm_square_bar_real_full["Ag"]
                + norm_test[3] * symm_comp_norm_square_bar_imag_full["Ag"]
            )

            norm_square_B1g = (
                norm_test[0] * symm_comp_norm_square_unbar_real_full["B1g"]
                + norm_test[1] * symm_comp_norm_square_unbar_imag_full["B1g"]
            )
            norm_square_B1g += (
                norm_test[2] * symm_comp_norm_square_bar_real_full["B1g"]
                + norm_test[3] * symm_comp_norm_square_bar_imag_full["B1g"]
            )

            norm_square_B2g = (
                norm_test[0] * symm_comp_norm_square_unbar_real_full["B2g"]
                + norm_test[1] * symm_comp_norm_square_unbar_imag_full["B2g"]
            )
            norm_square_B2g += (
                norm_test[2] * symm_comp_norm_square_bar_real_full["B2g"]
                + norm_test[3] * symm_comp_norm_square_bar_imag_full["B2g"]
            )

            norm_square_B3g = (
                norm_test[0] * symm_comp_norm_square_unbar_real_full["B3g"]
                + norm_test[1] * symm_comp_norm_square_unbar_imag_full["B3g"]
            )
            norm_square_B3g += (
                norm_test[2] * symm_comp_norm_square_bar_real_full["B3g"]
                + norm_test[3] * symm_comp_norm_square_bar_imag_full["B3g"]
            )

            # print("norm_square_Ag  = %20.16f" % (norm_square_Ag))
            # print("norm_square_B1g = %20.16f" % (norm_square_B1g))
            # print("norm_square_B2g = %20.16f" % (norm_square_B2g))
            # print("norm_square_B3g = %20.16f" % (norm_square_B3g))

            # c = 0.0
            # d = 0.0
            # e = 0.0
            # f = 0.0

            argmax_norm_test = numpy.argmax(entropy)
            if (
                norm_test[argmax_norm_test] < non_vanishing_tol
                or always_use_max_norm_comp
            ):
                argmax_norm_test = numpy.argmax(norm_test)
            # print(argmax_norm_test)

            dominant_irrep = None
            only_one_irrep = None

            if norm_square_Ag > 1e-14:
                if dominant_irrep is None and only_one_irrep is None:
                    dominant_irrep = "Ag"
                    only_one_irrep = True
                elif only_one_irrep and dominant_irrep != "Ag":
                    only_one_irrep = False
                    dominant_irrep = None
            if norm_square_B1g > 1e-14:
                if dominant_irrep is None and only_one_irrep is None:
                    dominant_irrep = "B1g"
                    only_one_irrep = True
                elif only_one_irrep and dominant_irrep != "B1g":
                    only_one_irrep = False
                    dominant_irrep = None
            if norm_square_B2g > 1e-14:
                if dominant_irrep is None and only_one_irrep is None:
                    dominant_irrep = "B2g"
                    only_one_irrep = True
                elif only_one_irrep and dominant_irrep != "B2g":
                    only_one_irrep = False
                    dominant_irrep = None
            if norm_square_B3g > 1e-14:
                if dominant_irrep is None and only_one_irrep is None:
                    dominant_irrep = "B3g"
                    only_one_irrep = True
                elif only_one_irrep and dominant_irrep != "B3g":
                    only_one_irrep = False
                    dominant_irrep = None

            if only_one_irrep:
                print("warning : only one irrep dominant irrep is ", dominant_irrep)

            if argmax_norm_test == 0:
                if norm_square_Ag > 1e-14:
                    c = numpy.sqrt(
                        norm_test[0]
                        * symm_comp_norm_square_unbar_real_full["Ag"]
                        / norm_square_Ag
                    )
                    if dominant_irrep is None and only_one_irrep is None:
                        dominant_irrep = "Ag"
                        only_one_irrep = True
                    elif only_one_irrep and dominant_irrep != "Ag":
                        only_one_irrep = False
                        dominant_irrep = None
                if norm_square_B1g > 1e-14:
                    d = -numpy.sqrt(
                        norm_test[0]
                        * symm_comp_norm_square_unbar_real_full["B1g"]
                        / norm_square_B1g
                    )
                if norm_square_B2g > 1e-14:
                    e = numpy.sqrt(
                        norm_test[0]
                        * symm_comp_norm_square_unbar_real_full["B2g"]
                        / norm_square_B2g
                    )
                if norm_square_B3g > 1e-14:
                    f = -numpy.sqrt(
                        norm_test[0]
                        * symm_comp_norm_square_unbar_real_full["B3g"]
                        / norm_square_B3g
                    )
            elif argmax_norm_test == 1:
                if norm_square_Ag > 1e-14:
                    d = numpy.sqrt(
                        norm_test[1]
                        * symm_comp_norm_square_unbar_imag_full["Ag"]
                        / norm_square_Ag
                    )
                if norm_square_B1g > 1e-14:
                    c = numpy.sqrt(
                        norm_test[1]
                        * symm_comp_norm_square_unbar_imag_full["B1g"]
                        / norm_square_B1g
                    )
                if norm_square_B2g > 1e-14:
                    f = numpy.sqrt(
                        norm_test[1]
                        * symm_comp_norm_square_unbar_imag_full["B2g"]
                        / norm_square_B2g
                    )
                if norm_square_B3g > 1e-14:
                    e = numpy.sqrt(
                        norm_test[1]
                        * symm_comp_norm_square_unbar_imag_full["B3g"]
                        / norm_square_B3g
                    )
            elif argmax_norm_test == 2:
                if norm_square_Ag > 1e-14:
                    e = -numpy.sqrt(
                        norm_test[2]
                        * symm_comp_norm_square_bar_real_full["Ag"]
                        / norm_square_Ag
                    )
                if norm_square_B1g > 1e-14:
                    f = -numpy.sqrt(
                        norm_test[2]
                        * symm_comp_norm_square_bar_real_full["B1g"]
                        / norm_square_B1g
                    )
                if norm_square_B2g > 1e-14:
                    c = numpy.sqrt(
                        norm_test[2]
                        * symm_comp_norm_square_bar_real_full["B2g"]
                        / norm_square_B2g
                    )
                if norm_square_B3g > 1e-14:
                    d = numpy.sqrt(
                        norm_test[2]
                        * symm_comp_norm_square_bar_real_full["B3g"]
                        / norm_square_B3g
                    )
            else:
                if norm_square_Ag > 1e-14:
                    f = numpy.sqrt(
                        norm_test[3]
                        * symm_comp_norm_square_bar_imag_full["Ag"]
                        / norm_square_Ag
                    )
                if norm_square_B1g > 1e-14:
                    e = -numpy.sqrt(
                        norm_test[3]
                        * symm_comp_norm_square_bar_imag_full["B1g"]
                        / norm_square_B1g
                    )
                if norm_square_B2g > 1e-14:
                    d = -numpy.sqrt(
                        norm_test[3]
                        * symm_comp_norm_square_bar_imag_full["B2g"]
                        / norm_square_B2g
                    )
                if norm_square_B3g > 1e-14:
                    c = numpy.sqrt(
                        norm_test[3]
                        * symm_comp_norm_square_bar_imag_full["B3g"]
                        / norm_square_B3g
                    )

            if (
                only_one_irrep
            ):  # 这种情况下 四种对称性只有 1 个有 不小于 0 的分量，此时 bar 或者 unbar 的 large comp 有一个数值上为 0 ! #
                use_real = False
                # if symm_comp_norm_square_unbar_real_full[dominant_irrep] > (1-pure_tol) and norm_test[0] > non_vanishing_tol:
                if norm_test[0] > non_vanishing_tol or norm_test[2] > non_vanishing_tol:
                    use_real = True
                print("use_real = ", use_real)

                if use_real:
                    symorb_unbar_irrep = symmorb_unbar_real_full[dominant_irrep]
                    symorb_bar_irrep = symmorb_bar_real_full[dominant_irrep]
                else:
                    symorb_unbar_irrep = symmorb_unbar_imag_full[dominant_irrep]
                    symorb_bar_irrep = symmorb_bar_imag_full[dominant_irrep]

                assert symorb_unbar_irrep is not None and symorb_bar_irrep is not None

                sum_mat = symorb_unbar_irrep - symorb_bar_irrep
                minus_mat = symorb_unbar_irrep + symorb_bar_irrep

                max_norm_sum = numpy.max(numpy.abs(sum_mat))
                max_norm_minus = numpy.max(numpy.abs(minus_mat))

                minus_sign = max_norm_sum > max_norm_minus

                print(
                    "symorb minus sign ? ",
                    minus_sign,
                    " with max_diff",
                    min(max_norm_sum, max_norm_minus),
                )

                # real rotation is enough #

                if use_real:
                    c_unbar = numpy.sqrt(
                        norm_test[0]
                        * symm_comp_norm_square_unbar_real_full[dominant_irrep]
                    )
                    c_bar = numpy.sqrt(
                        norm_test[2]
                        * symm_comp_norm_square_bar_real_full[dominant_irrep]
                    )
                else:
                    c_unbar = numpy.sqrt(
                        norm_test[1]
                        * symm_comp_norm_square_unbar_imag_full[dominant_irrep]
                    )
                    c_bar = numpy.sqrt(
                        norm_test[3]
                        * symm_comp_norm_square_bar_imag_full[dominant_irrep]
                    )

                if minus_sign:
                    c_bar = -c_bar

                norm_c = numpy.sqrt(c_unbar**2 + c_bar**2)
                c_unbar /= norm_c
                c_bar /= norm_c

                print("c_unbar = ", c_unbar, " c_bar = ", c_bar)

                c = 0.0
                d = 0.0
                e = 0.0
                f = 0.0

                if use_real:
                    if dominant_irrep == "Ag":
                        c = c_unbar
                        e = -c_bar
                    elif dominant_irrep == "B1g":
                        d = -c_unbar
                        f = -c_bar
                    elif dominant_irrep == "B2g":
                        e = c_unbar
                        c = c_bar
                    elif dominant_irrep == "B3g":
                        f = -c_unbar
                        d = c_bar
                    else:
                        raise ValueError("Impossible dominant irrep ", dominant_irrep)
                else:
                    if dominant_irrep == "Ag":
                        d = c_unbar
                        f = c_bar
                    elif dominant_irrep == "B1g":
                        c = c_unbar
                        e = -c_bar
                    elif dominant_irrep == "B2g":
                        f = c_unbar
                        d = -c_bar
                    elif dominant_irrep == "B3g":
                        e = c_unbar
                        c = c_bar
                    else:
                        raise ValueError("Impossible dominant irrep ", dominant_irrep)

        else:

            norm_square_Au = (
                norm_test[0] * symm_comp_norm_square_unbar_real_full["Au"]
                + norm_test[1] * symm_comp_norm_square_unbar_imag_full["Au"]
            )
            norm_square_Au += (
                norm_test[2] * symm_comp_norm_square_bar_real_full["Au"]
                + norm_test[3] * symm_comp_norm_square_bar_imag_full["Au"]
            )

            norm_square_B1u = (
                norm_test[0] * symm_comp_norm_square_unbar_real_full["B1u"]
                + norm_test[1] * symm_comp_norm_square_unbar_imag_full["B1u"]
            )
            norm_square_B1u += (
                norm_test[2] * symm_comp_norm_square_bar_real_full["B1u"]
                + norm_test[3] * symm_comp_norm_square_bar_imag_full["B1u"]
            )

            norm_square_B2u = (
                norm_test[0] * symm_comp_norm_square_unbar_real_full["B2u"]
                + norm_test[1] * symm_comp_norm_square_unbar_imag_full["B2u"]
            )
            norm_square_B2u += (
                norm_test[2] * symm_comp_norm_square_bar_real_full["B2u"]
                + norm_test[3] * symm_comp_norm_square_bar_imag_full["B2u"]
            )

            norm_square_B3u = (
                norm_test[0] * symm_comp_norm_square_unbar_real_full["B3u"]
                + norm_test[1] * symm_comp_norm_square_unbar_imag_full["B3u"]
            )
            norm_square_B3u += (
                norm_test[2] * symm_comp_norm_square_bar_real_full["B3u"]
                + norm_test[3] * symm_comp_norm_square_bar_imag_full["B3u"]
            )

            # print("norm_square_Ag  = %20.16f" % (norm_square_Ag))
            # print("norm_square_B1g = %20.16f" % (norm_square_B1g))
            # print("norm_square_B2g = %20.16f" % (norm_square_B2g))
            # print("norm_square_B3g = %20.16f" % (norm_square_B3g))
            # c = 0.0
            # d = 0.0
            # e = 0.0
            # f = 0.0

            argmax_norm_test = numpy.argmax(entropy)
            if (
                norm_test[argmax_norm_test] < non_vanishing_tol
                or always_use_max_norm_comp
            ):
                argmax_norm_test = numpy.argmax(norm_test)
            # argmax_norm_test = numpy.argmax(norm_test)
            # print(argmax_norm_test)

            dominant_irrep = None
            only_one_irrep = None

            if norm_square_Au > 1e-14:
                if dominant_irrep is None and only_one_irrep is None:
                    dominant_irrep = "Au"
                    only_one_irrep = True
                elif only_one_irrep and dominant_irrep != "Au":
                    only_one_irrep = False
                    dominant_irrep = None
            if norm_square_B1u > 1e-14:
                if dominant_irrep is None and only_one_irrep is None:
                    dominant_irrep = "B1u"
                    only_one_irrep = True
                elif only_one_irrep and dominant_irrep != "B1u":
                    only_one_irrep = False
                    dominant_irrep = None
            if norm_square_B2u > 1e-14:
                if dominant_irrep is None and only_one_irrep is None:
                    dominant_irrep = "B2u"
                    only_one_irrep = True
                elif only_one_irrep and dominant_irrep != "B2u":
                    only_one_irrep = False
                    dominant_irrep = None
            if norm_square_B3u > 1e-14:
                if dominant_irrep is None and only_one_irrep is None:
                    dominant_irrep = "B3u"
                    only_one_irrep = True
                elif only_one_irrep and dominant_irrep != "B3u":
                    only_one_irrep = False
                    dominant_irrep = None

            if only_one_irrep:
                print("warning : only one irrep dominant irrep is ", dominant_irrep)

            if argmax_norm_test == 0:
                if norm_square_Au > 1e-14:
                    c = numpy.sqrt(
                        norm_test[0]
                        * symm_comp_norm_square_unbar_real_full["Au"]
                        / norm_square_Au
                    )
                if norm_square_B1u > 1e-14:
                    d = -numpy.sqrt(
                        norm_test[0]
                        * symm_comp_norm_square_unbar_real_full["B1u"]
                        / norm_square_B1u
                    )
                if norm_square_B2u > 1e-14:
                    e = numpy.sqrt(
                        norm_test[0]
                        * symm_comp_norm_square_unbar_real_full["B2u"]
                        / norm_square_B2u
                    )
                if norm_square_B3u > 1e-14:
                    f = -numpy.sqrt(
                        norm_test[0]
                        * symm_comp_norm_square_unbar_real_full["B3u"]
                        / norm_square_B3u
                    )
            elif argmax_norm_test == 1:
                if norm_square_Au > 1e-14:
                    d = numpy.sqrt(
                        norm_test[1]
                        * symm_comp_norm_square_unbar_imag_full["Au"]
                        / norm_square_Au
                    )
                if norm_square_B1u > 1e-14:
                    c = numpy.sqrt(
                        norm_test[1]
                        * symm_comp_norm_square_unbar_imag_full["B1u"]
                        / norm_square_B1u
                    )
                if norm_square_B2u > 1e-14:
                    f = numpy.sqrt(
                        norm_test[1]
                        * symm_comp_norm_square_unbar_imag_full["B2u"]
                        / norm_square_B2u
                    )
                if norm_square_B3u > 1e-14:
                    e = numpy.sqrt(
                        norm_test[1]
                        * symm_comp_norm_square_unbar_imag_full["B3u"]
                        / norm_square_B3u
                    )
            elif argmax_norm_test == 2:
                if norm_square_Au > 1e-14:
                    e = -numpy.sqrt(
                        norm_test[2]
                        * symm_comp_norm_square_bar_real_full["Au"]
                        / norm_square_Au
                    )
                if norm_square_B1u > 1e-14:
                    f = -numpy.sqrt(
                        norm_test[2]
                        * symm_comp_norm_square_bar_real_full["B1u"]
                        / norm_square_B1u
                    )
                if norm_square_B2u > 1e-14:
                    c = numpy.sqrt(
                        norm_test[2]
                        * symm_comp_norm_square_bar_real_full["B2u"]
                        / norm_square_B2u
                    )
                if norm_square_B3u > 1e-14:
                    d = numpy.sqrt(
                        norm_test[2]
                        * symm_comp_norm_square_bar_real_full["B3u"]
                        / norm_square_B3u
                    )
            else:
                if norm_square_Au > 1e-14:
                    f = numpy.sqrt(
                        norm_test[3]
                        * symm_comp_norm_square_bar_imag_full["Au"]
                        / norm_square_Au
                    )
                if norm_square_B1u > 1e-14:
                    e = -numpy.sqrt(
                        norm_test[3]
                        * symm_comp_norm_square_bar_imag_full["B1u"]
                        / norm_square_B1u
                    )
                if norm_square_B2u > 1e-14:
                    d = -numpy.sqrt(
                        norm_test[3]
                        * symm_comp_norm_square_bar_imag_full["B2u"]
                        / norm_square_B2u
                    )
                if norm_square_B3u > 1e-14:
                    c = numpy.sqrt(
                        norm_test[3]
                        * symm_comp_norm_square_bar_imag_full["B3u"]
                        / norm_square_B3u
                    )

            if only_one_irrep:
                use_real = False
                # if symm_comp_norm_square_unbar_real_full[dominant_irrep] > (1-pure_tol) and norm_test[0] > non_vanishing_tol:
                if norm_test[0] > non_vanishing_tol or norm_test[2] > non_vanishing_tol:
                    use_real = True
                print("use_real = ", use_real)

                if use_real:
                    symorb_unbar_irrep = symmorb_unbar_real_full[dominant_irrep]
                    symorb_bar_irrep = symmorb_bar_real_full[dominant_irrep]
                else:
                    symorb_unbar_irrep = symmorb_unbar_imag_full[dominant_irrep]
                    symorb_bar_irrep = symmorb_bar_imag_full[dominant_irrep]

                assert symorb_unbar_irrep is not None and symorb_bar_irrep is not None

                sum_mat = symorb_unbar_irrep - symorb_bar_irrep
                minus_mat = symorb_unbar_irrep + symorb_bar_irrep

                max_norm_sum = numpy.max(numpy.abs(sum_mat))
                max_norm_minus = numpy.max(numpy.abs(minus_mat))

                minus_sign = max_norm_sum > max_norm_minus

                print(
                    "symorb minus sign ? ",
                    minus_sign,
                    " with max_diff",
                    min(max_norm_sum, max_norm_minus),
                )

                # real rotation is enough #

                if use_real:
                    c_unbar = numpy.sqrt(
                        norm_test[0]
                        * symm_comp_norm_square_unbar_real_full[dominant_irrep]
                    )
                    c_bar = numpy.sqrt(
                        norm_test[2]
                        * symm_comp_norm_square_bar_real_full[dominant_irrep]
                    )
                else:
                    c_unbar = numpy.sqrt(
                        norm_test[1]
                        * symm_comp_norm_square_unbar_imag_full[dominant_irrep]
                    )
                    c_bar = numpy.sqrt(
                        norm_test[3]
                        * symm_comp_norm_square_bar_imag_full[dominant_irrep]
                    )

                if minus_sign:
                    c_bar = -c_bar

                norm_c = numpy.sqrt(c_unbar**2 + c_bar**2)
                c_unbar /= norm_c
                c_bar /= norm_c

                print("c_unbar = ", c_unbar, " c_bar = ", c_bar)

                c = 0.0
                d = 0.0
                e = 0.0
                f = 0.0

                if use_real:
                    if dominant_irrep == "Au":
                        c = c_unbar
                        e = -c_bar
                    elif dominant_irrep == "B1u":
                        d = -c_unbar
                        f = -c_bar
                    elif dominant_irrep == "B2u":
                        e = c_unbar
                        c = c_bar
                    elif dominant_irrep == "B3u":
                        f = -c_unbar
                        d = c_bar
                    else:
                        raise ValueError("Impossible dominant irrep ", dominant_irrep)
                else:
                    if dominant_irrep == "Au":
                        d = c_unbar
                        f = c_bar
                    elif dominant_irrep == "B1u":
                        c = c_unbar
                        e = -c_bar
                    elif dominant_irrep == "B2u":
                        f = c_unbar
                        d = -c_bar
                    elif dominant_irrep == "B3u":
                        e = c_unbar
                        c = c_bar
                    else:
                        raise ValueError("Impossible dominant irrep ", dominant_irrep)

        # print(c, d, e, f)
        # print(c*c+d*d+e*e+f*f)

        norm_cdef = c * c + d * d + e * e + f * f

        if abs(norm_cdef - 1.0) > non_vanishing_tol:
            print("warning : norm_cdef = %20.16f not close to 1.0" % (norm_cdef))

        c /= numpy.sqrt(norm_cdef)
        d /= numpy.sqrt(norm_cdef)
        e /= numpy.sqrt(norm_cdef)
        f /= numpy.sqrt(norm_cdef)

        a = c + d * 1.0j
        b = e + f * 1.0j

        U_dag = numpy.zeros((2, 2), dtype=numpy.complex128)
        U_dag[0, 0] = numpy.conj(a)
        U_dag[1, 1] = a
        U_dag[0, 1] = numpy.conj(b)
        U_dag[1, 0] = -b

        # print(U_dag)

        mo_new_unbar = res[:, i] * numpy.conj(a) + res[:, i + 1] * (-b)
        mo_new_bar = res[:, i] * numpy.conj(b) + res[:, i + 1] * a
        res[:, i], res[:, i + 1] = mo_new_unbar, mo_new_bar

    return res


def purify_pes_d2h(mol, mf, pure_tol=1e-10, non_vanishing_tol=1e-8, max_iter=8):

    pes_mocoeff = mf.mo_coeff[:, mol.nao_2c() :].copy()
    adapted_spinor, _ = check_spinor_adapted(mol, pes_mocoeff)
    pes_mocoeff = purify_spinor_with_dominant_irrep(
        mol, pes_mocoeff, adapted_spinor, pure_tol, non_vanishing_tol
    )
    adapted_spinor, spinor_dominant_comp = check_spinor_adapted(mol, pes_mocoeff)

    ##############################
    # cache non-adapted spinor #
    ##############################

    spinor_cached = []
    spinor_cached_dominant_comp = []

    for _ in range(pes_mocoeff.shape[1] // 2):
        spinor_cached.append(None)
        spinor_cached_dominant_comp.append(1.0)
    for i in range(pes_mocoeff.shape[1] // 2):
        if not adapted_spinor[i]:
            spinor_cached[i] = pes_mocoeff[:, 2 * i : 2 * i + 2].copy()
            spinor_cached_dominant_comp[i] = spinor_dominant_comp[i]

    ##############################
    # do the work #
    ##############################

    for _ in range(max_iter):
        if not all(adapted_spinor):
            pes_mocoeff = purify_spinor_with_mixed_irrep(
                mol, pes_mocoeff, adapted_spinor, pure_tol, non_vanishing_tol
            )
            adapted_spinor, spinor_dominant_comp = check_spinor_adapted(
                mol, pes_mocoeff
            )

            for i in range(pes_mocoeff.shape[1] // 2):
                if not adapted_spinor[i] and (
                    spinor_dominant_comp[i] > spinor_cached_dominant_comp[i]
                ):
                    spinor_cached[i] = pes_mocoeff[:, 2 * i : 2 * i + 2].copy()
                    spinor_cached_dominant_comp[i] = spinor_dominant_comp[i]
                elif adapted_spinor[i]:
                    spinor_cached[i] = None
                    spinor_cached_dominant_comp[i] = 1.0

    for _ in range(max_iter):
        if not all(adapted_spinor):
            pes_mocoeff = purify_spinor_with_mixed_irrep(
                mol, pes_mocoeff, adapted_spinor, pure_tol, non_vanishing_tol, True
            )
            adapted_spinor, spinor_dominant_comp = check_spinor_adapted(
                mol, pes_mocoeff
            )

            for i in range(pes_mocoeff.shape[1] // 2):
                if not adapted_spinor[i] and (
                    spinor_dominant_comp[i] > spinor_cached_dominant_comp[i]
                ):
                    spinor_cached[i] = pes_mocoeff[:, 2 * i : 2 * i + 2].copy()
                    spinor_cached_dominant_comp[i] = spinor_dominant_comp[i]
                elif adapted_spinor[i]:
                    spinor_cached[i] = None
                    spinor_cached_dominant_comp[i] = 1.0

    ############################################################
    # if not success use the cached one with maximal dominant comp #
    ############################################################

    if not all(adapted_spinor):
        for i in range(pes_mocoeff.shape[1] // 2):
            if not adapted_spinor[i]:
                print(
                    "Use cached spinor %3d with dominant comp %20.16f"
                    % (i, spinor_cached_dominant_comp[i])
                )
                pes_mocoeff[:, 2 * i : 2 * i + 2] = spinor_cached[i].copy()

    check_spinor_adapted(mol, pes_mocoeff)

    return pes_mocoeff


import numpy as np


def warshall(A):
    """
    使用 Warshall 算法计算有向图的传递闭包（可达矩阵）。

    参数:
        A: 二维 0/1 矩阵，A[i][j] = 1 表示 i 到 j 有直接边
    返回:
        R: 布尔矩阵，R[i][j] = True 表示 i 到 j 有通路
    """
    A = np.asarray(A, dtype=bool)
    n = A.shape[0]

    R = A.copy()

    # 如果认为“点到自己总是可达的”，保留下面这行；
    # 如果只关心不同点之间，可以注释掉。
    np.fill_diagonal(R, True)

    for k in range(n):
        # 所有 i 能到 k，且 k 能到 j，则 i 能到 j
        R |= np.logical_and.outer(R[:, k], R[k, :])

    R = np.asarray(R, dtype=int)

    return R


def purify_pes_v2(mol, mf, vanish_tol=1e-10, align_tol=1e-5):  # much better ?

    print("------------ purify_pes_v2 -------------")

    pes_mocoeff = mf.mo_coeff[:, mol.nao_2c() :].copy()

    hcore_ao = mf.get_hcore()
    hcore_mo = mf.mo_coeff.conj().T @ hcore_ao @ mf.mo_coeff

    hcore_mo_imag = hcore_mo[mol.nao_2c() :, mol.nao_2c() :].imag
    if not numpy.allclose(hcore_mo_imag, 0.0, atol=1e-10):
        print(
            "warning : hcore_mo_imag not close to zero, max = ",
            numpy.max(numpy.abs(hcore_mo_imag)),
        )
        # print(
        #     "we cannot do the purification for d2h and its subgroup, since the hcore_mo is not real in pes space"
        # )
        # return pes_mocoeff
    hcore_mo = hcore_mo[mol.nao_2c() :, mol.nao_2c() :].real

    n2c = mol.nao_2c()

    # print 一轮 #

    aligned = numpy.zeros(
        (pes_mocoeff.shape[1] // 2, pes_mocoeff.shape[1] // 2), dtype=bool
    )

    for i in range(pes_mocoeff.shape[1] // 2):
        # for i in range(2):
        print("do i ", i)
        for j in range(i + 1, pes_mocoeff.shape[1] // 2):
            a1 = hcore_mo[2 * i, 2 * j]
            b1 = hcore_mo[2 * i, 2 * j + 1]
            a2 = hcore_mo[2 * i + 1, 2 * j + 1]
            b2 = hcore_mo[2 * i + 1, 2 * j]

            block = numpy.array([[a1, b1], [b2, a2]])
            # print("block %3d %3d" % (i, j), block)

            if abs(a1) > align_tol and abs(b1) < vanish_tol:
                aligned[i, j] = True
                aligned[j, i] = True
                # print("block %3d %3d aligned" % (i, j))
            else:
                if abs(a1) > align_tol:
                    if abs(b1) > vanish_tol:
                        if aligned[i, j]:
                            print(
                                "warning : block %3d %3d has non-vanishing off-diagonal elements % 12.6e"
                                % (i, j, abs(b1))
                            )
                        else:
                            print("align block %3d %3d" % (i, j))
                            ratio_ef = b1 / a1
                            e = 1.0 / numpy.sqrt(1 + abs(ratio_ef) ** 2)
                            f = ratio_ef * e
                            print(e, f)
                            mo_new = pes_mocoeff[:, 2 * j : 2 * j + 2] @ numpy.array(
                                [[e, -f], [f, e]]
                            )
                            pes_mocoeff[:, 2 * j : 2 * j + 2] = mo_new
                            aligned[j, i] = True
                            aligned[i, j] = True

        hcore_mo = pes_mocoeff.conj().T @ hcore_ao @ pes_mocoeff
        hcore_mo = hcore_mo.real
        aligned = warshall(aligned)

    aligned = warshall(aligned)

    for i in range(pes_mocoeff.shape[1] // 2):
        for j in range(i + 1, pes_mocoeff.shape[1] // 2):
            if aligned[i, j]:
                a1 = hcore_mo[2 * i, 2 * j]
                b1 = hcore_mo[2 * i, 2 * j + 1]
                a2 = hcore_mo[2 * i + 1, 2 * j + 1]
                b2 = hcore_mo[2 * i + 1, 2 * j]

                block = numpy.array([[a1, b1], [b2, a2]])
                # print("aligned block %3d %3d" % (i, j), block)
                # print("abs b1     = ", abs(b1))
                # print("vanish_tol = ", vanish_tol)
                if abs(b1) > vanish_tol:
                    print(
                        "warning : aligned block %3d %3d has non-vanishing off-diagonal elements %12.6e"
                        % (i, j, abs(b1))
                    )

    print("------------ purify_pes_v2 done -------------")

    return pes_mocoeff


# purify for d2h and its subgroup #

if __name__ == "__main__":

    from Lalpha_ao2spinor import get_ao2L_spinor
    from pyscf.symm import geom
    import integral_4C

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

    purify_pes_v2(mol, mf)

    hcore_ao = mf.get_hcore()
    hcore_mo = mf.mo_coeff.conj().T @ hcore_ao @ mf.mo_coeff
    hcore_mo_pes = hcore_mo[mol.nao_2c() :, mol.nao_2c() :]
    for i in range(hcore_mo_pes.shape[0]):
        for j in range(hcore_mo_pes.shape[1]):
            if abs(i - j) % 2 != 0:
                if abs(hcore_mo_pes[i, j]) > 1e-8:
                    print(
                        "warning : hcore_mo_pes[%3d,%3d] = %20.16e not zero"
                        % (i, j, hcore_mo_pes[i, j])
                    )
    # exit(1)

    res = purify_pes_d2h(mol, mf, non_vanishing_tol=1e-12, max_iter=8)
    n2c = mol.nao_2c()
    # mo_energy = mf.mo_energy[n2c:]
    mo_coeff = mf.mo_coeff.copy()
    mo_coeff[:, n2c:] = res
    mf.mo_coeff = mo_coeff
    # mf.mo_energy = mo_energy
    res = purify_pes_v2(mol, mf)

    check_spinor_adapted(mol, res)

    # exit(1)

    # 检查轨道能是否出问题 #

    n2c = mol.nao_2c()
    mo_energy = mf.mo_energy[n2c:]
    mo_coeff = mf.mo_coeff.copy()
    mo_coeff[:, n2c:] = res
    mf.mo_coeff = mo_coeff
    mf.mo_energy = mo_energy

    fock_ao = mf.get_fock()
    ovlp_ao = mf.get_ovlp()

    # print(fock_ao.shape)

    ovlp_mo = mo_coeff.conj().T @ ovlp_ao @ mo_coeff
    fock_mo = mo_coeff.conj().T @ fock_ao @ mo_coeff

    fock_mo_offdg = fock_mo - numpy.diag(numpy.diag(fock_mo))
    # assert numpy.allclose(fock_mo_offdg[n2c:, n2c:], 0.0)
    print("max abs offdg %20.16e" % (numpy.max(numpy.abs(fock_mo_offdg[n2c:, n2c:]))))
    diff_ovlp = ovlp_mo - numpy.eye(ovlp_mo.shape[0])
    print("max abs diff ovlp %20.16e" % (numpy.max(numpy.abs(diff_ovlp))))
    # exit(1)

    mo_energy2 = numpy.diag(fock_mo)[n2c:]

    for x, y in zip(mo_energy, mo_energy2):
        print("%24.16f %24.16f" % (x, y))

    hcore_ao = mf.get_hcore()
    hcore_mo = mf.mo_coeff.conj().T @ hcore_ao @ mf.mo_coeff
    hcore_mo_pes = hcore_mo[mol.nao_2c() :, mol.nao_2c() :]
    for i in range(hcore_mo_pes.shape[0]):
        for j in range(hcore_mo_pes.shape[1]):
            if abs(i - j) % 2 != 0:
                if abs(hcore_mo_pes[i, j]) > 1e-8:
                    print(
                        "warning : hcore_mo_pes[%3d,%3d] = %20.16e not zero"
                        % (i, j, hcore_mo_pes[i, j])
                    )

    # dump 4C 积分 #

    coulomb, breit = integral_4C.FCIDUMP_Rela4C(
        mol, mf, False, npes=28, filename="FCIDUMP_C2H4_bench", debug=True
    )

    assert breit is None
    print(coulomb.shape)

    coulomb_real = coulomb.real
    coulomb_imag = coulomb.imag

    print(numpy.allclose(coulomb_imag, 0.0))
    print("max imag %20.16e" % (numpy.max(numpy.abs(coulomb_imag))))

    integral_4C.FCIDUMP_Rela4C(
        mol, mf, False, npes=28, filename="FCIDUMP_C2H4", IsComplex=False
    )

    for i in range(28):
        for j in range(28):
            for k in range(28):
                for l in range(28):
                    n_odd = 0
                    if i % 2 == 1:
                        n_odd += 1
                    if j % 2 == 1:
                        n_odd += 1
                    if k % 2 == 1:
                        n_odd += 1
                    if l % 2 == 1:
                        n_odd += 1
                    if n_odd % 2 == 1:
                        if abs(coulomb_real[i, j, k, l]) >= 1e-7:
                            print(
                                "warning : coulomb_real[%3d,%3d,%3d,%3d] = %20.16e not zero"
                                % (i, j, k, l, coulomb_real[i, j, k, l])
                            )
