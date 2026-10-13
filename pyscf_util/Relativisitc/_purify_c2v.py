from pyscf import gto, scf, lib, symm
import numpy
from functools import reduce
from pyscf import lib
from pyscf import tools
import pyscf

from _purify_d2h import label_orb_symm


########################################
# 检查 spinor 是否 adapted #
########################################


def check_spinor_adapted(mol, pes_mocoeff, pure_tol=1e-10, non_vanishing_tol=1e-8):

    print(" ---------------------- check_spinor_adapted ---------------------- ")

    ALLOWED_SYM = [["A1"], ["A2"], ["B1"], ["B2"]]
    ALLOWED_SYM_BARRED = [
        ["B1"],
        ["B2"],
        ["A1"],
        ["A2"],
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

        for icomp in range(4):

            # print("icomp = ", icomp)

            norm_comp = norm_lst[icomp]
            norm_comp_bar = norm_lst_bar[icomp]

            # print(norm_comp, norm_comp_bar)

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

                # print("unbar ", symm_comp_norm_square)

                argmax_symm_comp = numpy.argmax(numpy.array(symm_comp_norm_square))
                argmax_symm_comp_value = symm_comp_norm_square[argmax_symm_comp]
                argmax_symm_comp_irrep_name = mol.irrep_name[argmax_symm_comp]

                if argmax_symm_comp_value < (1 - pure_tol):
                    print(
                        "comp %d unbarred is not pure since no single irrep dominates with largest comp %20.16f"
                        % (icomp, argmax_symm_comp_value)
                    )
                    is_all_adapted = False
                    result_dominant_comp[i // 2] = argmax_symm_comp_value

                if is_all_adapted:
                    if argmax_symm_comp_irrep_name not in ALLOWED_SYM[icomp]:
                        print(
                            "comp %d unbarred is not pure since the irrep %s is not in %s"
                            % (icomp, argmax_symm_comp_irrep_name, ALLOWED_SYM[icomp])
                        )
                        is_all_adapted = False
                        result_dominant_comp[i // 2] = 0.0

            if not is_all_adapted:
                break

            if bar_non_vanishing:
                symm_comp_norm_square_bar, _ = label_orb_symm(
                    mol,
                    mol.irrep_name,
                    mol.symm_orb,
                    mo_coeff_lst_bar[icomp].reshape(-1, 1),
                    ovlp,
                )

                # print("bar   ", symm_comp_norm_square_bar)

                argmax_symm_comp = numpy.argmax(numpy.array(symm_comp_norm_square_bar))
                argmax_symm_comp_value = symm_comp_norm_square_bar[argmax_symm_comp]
                argmax_symm_comp_irrep_name = mol.irrep_name[argmax_symm_comp]

                if argmax_symm_comp_value < (1 - pure_tol):
                    print(
                        "comp %d barred is not pure since no single irrep dominates with largest comp %20.16f"
                        % (icomp, argmax_symm_comp_value)
                    )
                    is_all_adapted = False
                    result_dominant_comp[i // 2] = argmax_symm_comp_value

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
                        is_all_adapted = False
                        result_dominant_comp[i // 2] = 0.0

            if not is_all_adapted:
                break

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

    ALLOWED_SYM = [["A1"], ["A2"], ["B1"], ["B2"]]
    ALLOWED_SYM_BARRED = [
        ["B1"],
        ["B2"],
        ["A1"],
        ["A2"],
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

    IRREP_NAME_C2V_FULL = ["A1", "A2", "B1", "B2"]

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

        # IsG = None
        # symorb = None
        symorb = {}
        for x in IRREP_NAME_C2V_FULL:
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
            for x in IRREP_NAME_C2V_FULL:
                if x not in symm_comp_norm_square_unbar_real_full:
                    symm_comp_norm_square_unbar_real_full[x] = 0.0
                    symmorb_unbar_real_full[x] = None

        else:
            symm_comp_norm_square_unbar_real_full = {}
            for x in IRREP_NAME_C2V_FULL:
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
            for x in IRREP_NAME_C2V_FULL:
                if x not in symm_comp_norm_square_unbar_imag_full:
                    symm_comp_norm_square_unbar_imag_full[x] = 0.0
                    symmorb_unbar_imag_full[x] = None

        else:

            symm_comp_norm_square_unbar_imag_full = {}
            for x in IRREP_NAME_C2V_FULL:
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
            for x in IRREP_NAME_C2V_FULL:
                if x not in symm_comp_norm_square_bar_real_full:
                    symm_comp_norm_square_bar_real_full[x] = 0.0
                    symmorb_bar_real_full[x] = None

        else:

            symm_comp_norm_square_bar_real_full = {}
            for x in IRREP_NAME_C2V_FULL:
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
            for x in IRREP_NAME_C2V_FULL:
                if x not in symm_comp_norm_square_bar_imag_full:
                    symm_comp_norm_square_bar_imag_full[x] = 0.0
                    symmorb_bar_imag_full[x] = None

        else:

            symm_comp_norm_square_bar_imag_full = {}
            for x in IRREP_NAME_C2V_FULL:
                symm_comp_norm_square_bar_imag_full[x] = 0.0
            symmorb_bar_imag_full = None

        for name in IRREP_NAME_C2V_FULL:
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

            for x in IRREP_NAME_C2V_FULL:
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
            for x in IRREP_NAME_C2V_FULL:
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
            for x in IRREP_NAME_C2V_FULL:
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
            for x in IRREP_NAME_C2V_FULL:
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

        for x in IRREP_NAME_C2V_FULL:
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

        # if IsG:
        if True:

            norm_square_A1 = (
                norm_test[0] * symm_comp_norm_square_unbar_real_full["A1"]
                + norm_test[1] * symm_comp_norm_square_unbar_imag_full["A1"]
            )
            norm_square_A1 += (
                norm_test[2] * symm_comp_norm_square_bar_real_full["A1"]
                + norm_test[3] * symm_comp_norm_square_bar_imag_full["A1"]
            )

            norm_square_A2 = (
                norm_test[0] * symm_comp_norm_square_unbar_real_full["A2"]
                + norm_test[1] * symm_comp_norm_square_unbar_imag_full["A2"]
            )
            norm_square_A2 += (
                norm_test[2] * symm_comp_norm_square_bar_real_full["A2"]
                + norm_test[3] * symm_comp_norm_square_bar_imag_full["A2"]
            )

            norm_square_B1 = (
                norm_test[0] * symm_comp_norm_square_unbar_real_full["B1"]
                + norm_test[1] * symm_comp_norm_square_unbar_imag_full["B1"]
            )
            norm_square_B1 += (
                norm_test[2] * symm_comp_norm_square_bar_real_full["B1"]
                + norm_test[3] * symm_comp_norm_square_bar_imag_full["B1"]
            )

            norm_square_B2 = (
                norm_test[0] * symm_comp_norm_square_unbar_real_full["B2"]
                + norm_test[1] * symm_comp_norm_square_unbar_imag_full["B2"]
            )
            norm_square_B2 += (
                norm_test[2] * symm_comp_norm_square_bar_real_full["B2"]
                + norm_test[3] * symm_comp_norm_square_bar_imag_full["B2"]
            )

            # print("norm_square_A1 = %20.16f" % (norm_square_A1))
            # print("norm_square_A2 = %20.16f" % (norm_square_A2))
            # print("norm_square_B1 = %20.16f" % (norm_square_B1))
            # print("norm_square_B2 = %20.16f" % (norm_square_B2))

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

            if norm_square_A1 > 1e-14:
                if dominant_irrep is None and only_one_irrep is None:
                    dominant_irrep = "A1"
                    only_one_irrep = True
                elif only_one_irrep and dominant_irrep != "A1":
                    only_one_irrep = False
                    dominant_irrep = None
            if norm_square_A2 > 1e-14:
                if dominant_irrep is None and only_one_irrep is None:
                    dominant_irrep = "A2"
                    only_one_irrep = True
                elif only_one_irrep and dominant_irrep != "A2":
                    only_one_irrep = False
                    dominant_irrep = None
            if norm_square_B1 > 1e-14:
                if dominant_irrep is None and only_one_irrep is None:
                    dominant_irrep = "B1"
                    only_one_irrep = True
                elif only_one_irrep and dominant_irrep != "B1":
                    only_one_irrep = False
                    dominant_irrep = None
            if norm_square_B2 > 1e-14:
                if dominant_irrep is None and only_one_irrep is None:
                    dominant_irrep = "B2"
                    only_one_irrep = True
                elif only_one_irrep and dominant_irrep != "B2":
                    only_one_irrep = False
                    dominant_irrep = None

            if only_one_irrep:
                print("warning : only one irrep dominant irrep is ", dominant_irrep)

            if argmax_norm_test == 0:
                if norm_square_A1 > 1e-14:
                    c = numpy.sqrt(
                        norm_test[0]
                        * symm_comp_norm_square_unbar_real_full["A1"]
                        / norm_square_A1
                    )
                    if dominant_irrep is None and only_one_irrep is None:
                        dominant_irrep = "A1"
                        only_one_irrep = True
                    elif only_one_irrep and dominant_irrep != "A1":
                        only_one_irrep = False
                        dominant_irrep = None
                if norm_square_A2 > 1e-14:
                    d = -numpy.sqrt(
                        norm_test[0]
                        * symm_comp_norm_square_unbar_real_full["A2"]
                        / norm_square_A2
                    )
                if norm_square_B1 > 1e-14:
                    e = numpy.sqrt(
                        norm_test[0]
                        * symm_comp_norm_square_unbar_real_full["B1"]
                        / norm_square_B1
                    )
                if norm_square_B2 > 1e-14:
                    f = -numpy.sqrt(
                        norm_test[0]
                        * symm_comp_norm_square_unbar_real_full["B2"]
                        / norm_square_B2
                    )
            elif argmax_norm_test == 1:
                if norm_square_A1 > 1e-14:
                    d = numpy.sqrt(
                        norm_test[1]
                        * symm_comp_norm_square_unbar_imag_full["A1"]
                        / norm_square_A1
                    )
                if norm_square_A2 > 1e-14:
                    c = numpy.sqrt(
                        norm_test[1]
                        * symm_comp_norm_square_unbar_imag_full["A2"]
                        / norm_square_A2
                    )
                if norm_square_B1 > 1e-14:
                    f = numpy.sqrt(
                        norm_test[1]
                        * symm_comp_norm_square_unbar_imag_full["B1"]
                        / norm_square_B1
                    )
                if norm_square_B2 > 1e-14:
                    e = numpy.sqrt(
                        norm_test[1]
                        * symm_comp_norm_square_unbar_imag_full["B2"]
                        / norm_square_B2
                    )
            elif argmax_norm_test == 2:
                if norm_square_A1 > 1e-14:
                    e = -numpy.sqrt(
                        norm_test[2]
                        * symm_comp_norm_square_bar_real_full["A1"]
                        / norm_square_A1
                    )
                if norm_square_A2 > 1e-14:
                    f = -numpy.sqrt(
                        norm_test[2]
                        * symm_comp_norm_square_bar_real_full["A2"]
                        / norm_square_A2
                    )
                if norm_square_B1 > 1e-14:
                    c = numpy.sqrt(
                        norm_test[2]
                        * symm_comp_norm_square_bar_real_full["B1"]
                        / norm_square_B1
                    )
                if norm_square_B2 > 1e-14:
                    d = numpy.sqrt(
                        norm_test[2]
                        * symm_comp_norm_square_bar_real_full["B2"]
                        / norm_square_B2
                    )
            else:
                if norm_square_A1 > 1e-14:
                    f = numpy.sqrt(
                        norm_test[3]
                        * symm_comp_norm_square_bar_imag_full["A1"]
                        / norm_square_A1
                    )
                if norm_square_A2 > 1e-14:
                    e = -numpy.sqrt(
                        norm_test[3]
                        * symm_comp_norm_square_bar_imag_full["A2"]
                        / norm_square_A2
                    )
                if norm_square_B1 > 1e-14:
                    d = -numpy.sqrt(
                        norm_test[3]
                        * symm_comp_norm_square_bar_imag_full["B1"]
                        / norm_square_B1
                    )
                if norm_square_B2 > 1e-14:
                    c = numpy.sqrt(
                        norm_test[3]
                        * symm_comp_norm_square_bar_imag_full["B2"]
                        / norm_square_B2
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
                    if dominant_irrep == "A1":
                        c = c_unbar
                        e = -c_bar
                    elif dominant_irrep == "A2":
                        d = -c_unbar
                        f = -c_bar
                    elif dominant_irrep == "B1":
                        e = c_unbar
                        c = c_bar
                    elif dominant_irrep == "B2":
                        f = -c_unbar
                        d = c_bar
                    else:
                        raise ValueError("Impossible dominant irrep ", dominant_irrep)
                else:
                    if dominant_irrep == "A1":
                        d = c_unbar
                        f = c_bar
                    elif dominant_irrep == "A2":
                        c = c_unbar
                        e = -c_bar
                    elif dominant_irrep == "B1":
                        f = c_unbar
                        d = -c_bar
                    elif dominant_irrep == "B2":
                        e = c_unbar
                        c = c_bar
                    else:
                        raise ValueError("Impossible dominant irrep ", dominant_irrep)

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


def purify_pes_c2v(mol, mf, pure_tol=1e-10, non_vanishing_tol=1e-8, max_iter=8):

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
            spinor_cached[i] = pes_mocoeff[:, 2 * i : 2 * i + 2]
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
                    spinor_cached[i] = pes_mocoeff[:, 2 * i : 2 * i + 2]
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
                    spinor_cached[i] = pes_mocoeff[:, 2 * i : 2 * i + 2]
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
                pes_mocoeff[:, 2 * i : 2 * i + 2] = spinor_cached[i]

    return pes_mocoeff


# purify for d2h and its subgroup #

if __name__ == "__main__":

    from Lalpha_ao2spinor import get_ao2L_spinor
    from pyscf.symm import geom

    mol = gto.M(
        atom="""
O    0.000000   0.000000   0.000000
H    0.758000   0.000000   0.588000
H   -0.758000  -0.000000   0.588000
""",
        basis="unc-cc-pvdz-dk",
        # basis="3-21g",
        # basis="cc-pvdz-dk",
        verbose=5,
        charge=0,
        spin=0,
        symmetry="c2v",
        # unit='bohr'
    )

    # print(get_ao2L_spinor(mol))

    ao2Lalpha = get_ao2L_spinor(mol)
    ao2Lbeta = get_ao2L_spinor(mol, False)

    # print(mol._symm_)

    gpname, origin, axes = geom.detect_symm(mol._atom)
    print(gpname, origin, axes)
    # exit(1)

    # do scf #

    mol.build()
    mf = scf.dhf.RDHF(mol)
    mf.conv_tol = 1e-12
    mf.kernel()

    res = purify_pes_c2v(mol, mf, max_iter=8)

    # 检查轨道能是否出问题 #

    n2c = mol.nao_2c()
    mo_energy = mf.mo_energy[n2c:]
    mo_coeff = mf.mo_coeff.copy()
    mo_coeff[:, n2c:] = res

    fock_ao = mf.get_fock()
    print(fock_ao.shape)

    fock_mo = mo_coeff.conj().T @ fock_ao @ mo_coeff

    fock_mo_offdg = fock_mo - numpy.diag(numpy.diag(fock_mo))
    # assert numpy.allclose(fock_mo_offdg[n2c:, n2c:], 0.0)
    print("max abs offdg %20.16e" % (numpy.max(numpy.abs(fock_mo_offdg[n2c:, n2c:]))))

    mo_energy2 = numpy.diag(fock_mo)[n2c:]

    for x, y in zip(mo_energy, mo_energy2):
        print("%24.16f %24.16f" % (x, y))
