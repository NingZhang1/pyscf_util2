from pyscf_util.Analyzer.iCIPT2.analyzer_mrpt2 import *
from tabulate import tabulate

from openpyxl import Workbook


# 返回MRPT2总能量及其不确定度
def analysis_mrpt2_excel(
    ENPT2_RES,
    MRENPT2_RES,
    NEVPT2_RES,
    TASK,
    CMIN,
    wb,
    wb_title,
    print_pic=False,
    print_table=True,
    NPT=5,
    NLAST_REMOVE_LARGE=None,
    NLAST_REMOVE_SMALL=None,
    CMIN_MRENPT2=None,
    NPT_MRENTP2=None,
    only_extra_res=False,
):

    ws = wb.create_sheet(wb_title)
    # ws.title = wb_title

    # some global data #

    subspace_order = [(0, 1), (0, 2), (1, 0), (1, 1), (2, 0)]
    subspace_order_2_key = {
        (0, 1): "r",
        (0, 2): "rs",
        (1, 0): "i",
        (1, 1): "ir",
        (2, 0): "ij",
    }
    subspace_order2 = ["ijr", "rsi", "ijrs"]
    header = ["Cmin", "ept2", "r", "rs", "i", "ir", "ij", "ijr", "rsi", "ijrs"]
    subspace_order_print = ["r", "rs", "i", "ir", "ij", "ijr", "rsi", "ijrs"]

    if CMIN_MRENPT2 is None:
        CMIN_MRENPT2 = CMIN
    if NPT_MRENTP2 is None:
        NPT_MRENTP2 = NPT

    for i in range(len(CMIN_MRENPT2)):
        if abs(CMIN_MRENPT2[i] - CMIN[i]) > 1e-10:
            exit(1)
    len_cmin_mrenpt2 = len(CMIN_MRENPT2)

    return_res = {}

    for mole in TASK:

        return_res[mole] = {}
        ept_tot = 0.0
        ept_err = 0.0

        data_print = []
        DATA_EXTRA = {
            "i": {
                "ept": [],
                "etot": [],
            },
            "ij": {
                "ept": [],
                "etot": [],
            },
            "ir": {
                "ept": [],
                "etot": [],
            },
            "r": {
                "ept": [],
                "etot": [],
            },
            "rs": {
                "ept": [],
                "etot": [],
            },
            "ijr": {
                "ept": [],
                "etot": [],
            },
            "rsi": {
                "ept": [],
                "etot": [],
            },
            "ijrs": {
                "ept": [],
                "etot": [],
            },
        }

        # for cmin in CMIN:
        for cmin in CMIN_MRENPT2:

            # table res #

            data = [cmin]
            ept2 = ENPT2_RES[mole][cmin]["perturbation"]
            data.append(ept2)

            for key in subspace_order:
                e2 = MRENPT2_RES[(mole, cmin)][key]
                data.append(e2)

                key2 = subspace_order_2_key[key]
                DATA_EXTRA[key2]["ept"].append(ept2)
                DATA_EXTRA[key2]["etot"].append(e2)

            for key in subspace_order2:
                e2 = NEVPT2_RES[(mole)]["pc-NEVPT2"][cmin][key]["e"]
                data.append(e2)

                DATA_EXTRA[key]["ept"].append(ept2)
                DATA_EXTRA[key]["etot"].append(e2)

            data_print.append(data)

        for cmin in CMIN[len_cmin_mrenpt2:]:

            data = [cmin]
            ept2 = ENPT2_RES[mole][cmin]["perturbation"]
            data.append(ept2)

            for key in subspace_order:
                # e2 = MRENPT2_RES[(mole, cmin)][key]
                data.append(0.0)

                # key2 = subspace_order_2_key[key]
                # DATA_EXTRA[key2]["ept"].append(ept2)
                # DATA_EXTRA[key2]["etot"].append(e2)

            for key in subspace_order2:
                e2 = NEVPT2_RES[(mole)]["pc-NEVPT2"][cmin][key]["e"]
                data.append(e2)

                DATA_EXTRA[key]["ept"].append(ept2)
                DATA_EXTRA[key]["etot"].append(e2)

            data_print.append(data)

        # print(DATA_EXTRA)

        DATA_PRINT2 = {}
        DATA_PRINT3 = {}

        data_large_linear = ["large-linear", "energy"]
        data_large_linear_error = ["large-linear", "error"]
        data_large_quadratic = ["large-quadratic", "energy"]
        data_large_quadratic_error = ["large-quadratic", "error"]
        data_large_pade = ["large-pade", "energy"]
        data_large_pade_error = ["large-pade", "error"]

        data_small_linear = ["small-linear", "energy"]
        data_small_linear_error = ["small-linear", "error"]
        data_small_quadratic = ["small-quadratic", "energy"]
        data_small_quadratic_error = ["small-quadratic", "error"]
        data_small_pade = ["small-pade", "energy"]
        data_small_pade_error = ["small-pade", "error"]

        for key in subspace_order_print:
            if key in ["r", "rs", "i", "ir", "ij"]:
                DATA_PRINT2[key] = get_extra_res(
                    DATA_EXTRA[key]["ept"],
                    DATA_EXTRA[key]["etot"],
                    NPT_MRENTP2,
                    NLAST_REMOVE_SMALL,
                )
                DATA_PRINT3[key] = get_extra_res(
                    DATA_EXTRA[key]["ept"],
                    DATA_EXTRA[key]["etot"],
                    NPT_MRENTP2,
                    NLAST_REMOVE_LARGE,
                )
            else:
                DATA_PRINT2[key] = get_extra_res(
                    DATA_EXTRA[key]["ept"],
                    DATA_EXTRA[key]["etot"],
                    NPT,
                    NLAST_REMOVE_SMALL,
                )
                DATA_PRINT3[key] = get_extra_res(
                    DATA_EXTRA[key]["ept"],
                    DATA_EXTRA[key]["etot"],
                    NPT,
                    NLAST_REMOVE_LARGE,
                )

            data_large_linear.append(DATA_PRINT3[key]["weighted_linear_extra"])
            data_large_linear_error.append(DATA_PRINT3[key]["weighted_linear_error"])
            data_large_quadratic.append(DATA_PRINT3[key]["weighted_quadratic_extra"])
            data_large_quadratic_error.append(
                DATA_PRINT3[key]["weighted_quadratic_error"]
            )

            data_small_linear.append(DATA_PRINT2[key]["weighted_linear_extra"])
            data_small_linear_error.append(DATA_PRINT2[key]["weighted_linear_error"])
            data_small_quadratic.append(DATA_PRINT2[key]["weighted_quadratic_extra"])
            data_small_quadratic_error.append(
                DATA_PRINT2[key]["weighted_quadratic_error"]
            )

            data_large_pade.append(DATA_PRINT3[key]["weighted_pade_extra"])
            data_large_pade_error.append(DATA_PRINT3[key]["weighted_pade_error"])
            data_small_pade.append(DATA_PRINT2[key]["weighted_pade_extra"])
            data_small_pade_error.append(DATA_PRINT2[key]["weighted_pade_error"])

            ept_tot += DATA_PRINT2[key]["weighted_linear_extra"]
            ept_err += DATA_PRINT2[key]["weighted_linear_error"] ** 2

        data_print.append(data_large_linear)
        data_print.append(data_large_quadratic)
        data_print.append(data_large_pade)
        data_print.append(data_small_linear)
        data_print.append(data_small_quadratic)
        data_print.append(data_small_pade)

        data_print.append(data_large_linear_error)
        data_print.append(data_large_quadratic_error)
        data_print.append(data_large_pade_error)
        data_print.append(data_small_linear_error)
        data_print.append(data_small_quadratic_error)
        data_print.append(data_small_pade_error)

        if not only_extra_res:
            ws.append(header)
            for x in data_print:
                ws.append(x)
        else:
            ws.append(data_small_linear)

        # update return res #

        ept_err = np.sqrt(ept_err)
        return_res[mole] = {
            "ept": ept_tot,
            "err": ept_err,
        }

        # print #

        if print_pic or print_table:
            print(mole)

        if print_table:
            print(
                tabulate(data_print, headers=header, tablefmt="grid", floatfmt="15.8f")
            )

        if print_pic:
            print("extra with small cmin")
            draw_extra_pic(DATA_PRINT2, 2, 4, subspace_order_print, 24, 9)
            draw_extra_pic(
                DATA_PRINT2, 2, 4, subspace_order_print, 24, 9, use_quadratic=True
            )
            draw_extra_pic(
                DATA_PRINT2, 2, 4, subspace_order_print, 24, 9, use_pade=True
            )

            print("extra with large cmin")
            draw_extra_pic(DATA_PRINT3, 2, 4, subspace_order_print, 24, 9)
            draw_extra_pic(
                DATA_PRINT3, 2, 4, subspace_order_print, 24, 9, use_quadratic=True
            )
            draw_extra_pic(
                DATA_PRINT3, 2, 4, subspace_order_print, 24, 9, use_pade=True
            )
    
    return return_res


def analysis_mrpt2_2_excel(
    ENPT2_RES,
    MRENPT2_RES,
    NEVPT2_RES,
    TASK,
    CMIN,
    wb,
    wb_title,
    print_pic=False,
    print_table=True,
    NPT=5,
    NLAST_REMOVE_LARGE=None,
    NLAST_REMOVE_SMALL=None,
    CMIN_MRENPT2=None,
    NPT_MRENTP2=None,
):

    ws = wb.create_sheet(wb_title)
    # ws.title = wb_title

    # some global data #

    # subspace_order = [(0, 1), (0, 2), (1, 0), (1, 1), (2, 0)]
    # subspace_order_2_key = {
    #     (0, 1): "r",
    #     (0, 2): "rs",
    #     (1, 0): "i",
    #     (1, 1): "ir",
    #     (2, 0): "ij",
    # }
    subspace_order2 = ["ijr", "rsi", "ijrs"]
    header = ["Cmin", "ept2", "space", "ijr", "rsi", "ijrs"]
    subspace_order_print = ["space", "ijr", "rsi", "ijrs"]

    if CMIN_MRENPT2 is None:
        CMIN_MRENPT2 = CMIN
    if NPT_MRENTP2 is None:
        NPT_MRENTP2 = NPT

    for i in range(len(CMIN_MRENPT2)):
        if abs(CMIN_MRENPT2[i] - CMIN[i]) > 1e-10:
            exit(1)
    len_cmin_mrenpt2 = len(CMIN_MRENPT2)

    for mole in TASK:
        data_print = []
        DATA_EXTRA = {
            "space": {
                "ept": [],
                "etot": [],
            },
            "ijr": {
                "ept": [],
                "etot": [],
            },
            "rsi": {
                "ept": [],
                "etot": [],
            },
            "ijrs": {
                "ept": [],
                "etot": [],
            },
        }
        # for cmin in CMIN:

        for cmin in CMIN_MRENPT2:

            # table res #

            data = [cmin]
            ept2 = ENPT2_RES[mole][cmin]["perturbation"]
            data.append(ept2)

            # for key in subspace_order:
            e2 = MRENPT2_RES[(mole, cmin)]
            data.append(e2)

            # key2 = subspace_order_2_key[key]
            DATA_EXTRA["space"]["ept"].append(ept2)
            DATA_EXTRA["space"]["etot"].append(e2)

            for key in subspace_order2:
                e2 = NEVPT2_RES[(mole)]["pc-NEVPT2"][cmin][key]["e"]
                data.append(e2)

                DATA_EXTRA[key]["ept"].append(ept2)
                DATA_EXTRA[key]["etot"].append(e2)

            data_print.append(data)

        for cmin in CMIN[len_cmin_mrenpt2:]:

            data = [cmin]
            ept2 = ENPT2_RES[mole][cmin]["perturbation"]
            data.append(ept2)

            data.append(0.0)

            for key in subspace_order2:
                e2 = NEVPT2_RES[(mole)]["pc-NEVPT2"][cmin][key]["e"]
                data.append(e2)

                DATA_EXTRA[key]["ept"].append(ept2)
                DATA_EXTRA[key]["etot"].append(e2)

            data_print.append(data)

        # print(DATA_EXTRA)

        DATA_PRINT2 = {}
        DATA_PRINT3 = {}

        data_large_linear = ["large-linear", "energy"]
        data_large_linear_error = ["large-linear", "error"]
        data_large_quadratic = ["large-quadratic", "energy"]
        data_large_quadratic_error = ["large-quadratic", "error"]
        data_large_pade = ["large-pade", "energy"]
        data_large_pade_error = ["large-pade", "error"]

        data_small_linear = ["small-linear", "energy"]
        data_small_linear_error = ["small-linear", "error"]
        data_small_quadratic = ["small-quadratic", "energy"]
        data_small_quadratic_error = ["small-quadratic", "error"]
        data_small_pade = ["small-pade", "energy"]
        data_small_pade_error = ["small-pade", "error"]

        for key in subspace_order_print:
            if key == "space":
                DATA_PRINT2[key] = get_extra_res(
                    DATA_EXTRA[key]["ept"],
                    DATA_EXTRA[key]["etot"],
                    NPT_MRENTP2,
                    NLAST_REMOVE_SMALL,
                )
                DATA_PRINT3[key] = get_extra_res(
                    DATA_EXTRA[key]["ept"],
                    DATA_EXTRA[key]["etot"],
                    NPT_MRENTP2,
                    NLAST_REMOVE_LARGE,
                )
            else:
                DATA_PRINT2[key] = get_extra_res(
                    DATA_EXTRA[key]["ept"],
                    DATA_EXTRA[key]["etot"],
                    NPT,
                    NLAST_REMOVE_SMALL,
                )
                DATA_PRINT3[key] = get_extra_res(
                    DATA_EXTRA[key]["ept"],
                    DATA_EXTRA[key]["etot"],
                    NPT,
                    NLAST_REMOVE_LARGE,
                )

            data_large_linear.append(DATA_PRINT3[key]["weighted_linear_extra"])
            data_large_linear_error.append(DATA_PRINT3[key]["weighted_linear_error"])
            data_large_quadratic.append(DATA_PRINT3[key]["weighted_quadratic_extra"])
            data_large_quadratic_error.append(
                DATA_PRINT3[key]["weighted_quadratic_error"]
            )

            data_small_linear.append(DATA_PRINT2[key]["weighted_linear_extra"])
            data_small_linear_error.append(DATA_PRINT2[key]["weighted_linear_error"])
            data_small_quadratic.append(DATA_PRINT2[key]["weighted_quadratic_extra"])
            data_small_quadratic_error.append(
                DATA_PRINT2[key]["weighted_quadratic_error"]
            )

            data_large_pade.append(DATA_PRINT3[key]["weighted_pade_extra"])
            data_large_pade_error.append(DATA_PRINT3[key]["weighted_pade_error"])
            data_small_pade.append(DATA_PRINT2[key]["weighted_pade_extra"])
            data_small_pade_error.append(DATA_PRINT2[key]["weighted_pade_error"])

        data_print.append(data_large_linear)
        data_print.append(data_large_quadratic)
        data_print.append(data_large_pade)
        data_print.append(data_small_linear)
        data_print.append(data_small_quadratic)
        data_print.append(data_small_pade)

        data_print.append(data_large_linear_error)
        data_print.append(data_large_quadratic_error)
        data_print.append(data_large_pade_error)
        data_print.append(data_small_linear_error)
        data_print.append(data_small_quadratic_error)
        data_print.append(data_small_pade_error)

        ws.append(header)
        for x in data_print:
            # print(x)
            ws.append(x)

        # print #

        if print_pic or print_table:
            print(mole)

        if print_table:
            print(
                tabulate(data_print, headers=header, tablefmt="grid", floatfmt="15.8f")
            )

        if print_pic:
            print("extra with small cmin")
            draw_extra_pic(DATA_PRINT2, 2, 2, subspace_order_print, 12, 9)
            draw_extra_pic(
                DATA_PRINT2, 2, 2, subspace_order_print, 12, 9, use_quadratic=True
            )
            draw_extra_pic(
                DATA_PRINT2, 2, 2, subspace_order_print, 12, 9, use_pade=True
            )

            print("extra with large cmin")
            draw_extra_pic(DATA_PRINT3, 2, 2, subspace_order_print, 12, 9)
            draw_extra_pic(
                DATA_PRINT3, 2, 2, subspace_order_print, 12, 9, use_quadratic=True
            )
            draw_extra_pic(
                DATA_PRINT3, 2, 2, subspace_order_print, 12, 9, use_pade=True
            )
