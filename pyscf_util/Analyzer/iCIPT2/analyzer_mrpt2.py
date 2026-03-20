from pyscf_util.Analyzer.iCIPT2.analyzer import *
from tabulate import tabulate

############################################
## extract mrpt2/nevpt2 ##
############################################


def extract_mrenpt2_res(filename: str):
    with open(filename, "r") as file:
        content = file.read()
        if "--------------------- MRPT2 Driver End ---------------------" in content:
            matrix_elements = {}
            for line in content.splitlines():
                if re.match(r"\(\s*\d+,\s*\d+\s*\)\s*\|\s*-?\d+\.\d+", line):
                    indices, value = line.split("|")
                    i, j = map(int, re.findall(r"\d+", indices))
                    value = float(value.strip())
                    if not (
                        (i == 1 and j == 2)
                        or (i == 2 and j == 1)
                        or (i == 2 and j == 2)
                    ):
                        matrix_elements[(i, j)] = value
            return matrix_elements
        else:
            print(f"Error: File {filename} does not contain the required string.")


############################################
## extract nevpt2s ##
############################################

############################################
## case I one qmin extract both ept and etot
############################################


def extract_nevpt2s_old_type(filename: str):
    order = [
        (0, 1),
        (0, 2),
        (1, 0),
        (1, 1),
        (2, 0),
    ]
    with open(filename, "r") as file:
        content = file.read()
        if "--------------------- MRPT2 Driver End ---------------------" in content:
            pattern = r"perturbation.*\n\s*[^/]*/([^/]*)/"
            match = re.findall(pattern, content)
            assert len(match) == len(order)
            res = {}
            for i in range(len(order)):
                res[order[i]] = {"ept": float(match[i]), "etot": 0.0}

            # find etot #

            for line in content.splitlines():
                if re.match(r"\(\s*\d+,\s*\d+\s*\)\s*\|\s*-?\d+\.\d+", line):
                    indices, value = line.split("|")
                    i, j = map(int, re.findall(r"\d+", indices))
                    value = float(value.strip())
                    if not (
                        (i == 1 and j == 2)
                        or (i == 2 and j == 1)
                        or (i == 2 and j == 2)
                    ):
                        res[(i, j)]["etot"] = value

            return res

        else:
            print(f"Error: File {filename} does not contain the required string.")


def extract_nevpt2s_new_type(filename: str, with_LCUA=False, with_full_gFock=False):

    raw_data_factor = 10
    if with_full_gFock:
        with_LCUA = True
    if with_full_gFock:
        raw_data_factor += 5
    if with_LCUA:
        raw_data_factor += 5

    with open(filename, "r") as file:
        content = file.read()
        if "--------------------- MRPT2 Driver End ---------------------" in content:

            res_tmp = []

            pattern = r"Qmin\s*=\s*((?:\d+\.\d+e-?\d+\s*)+)"
            matches = re.findall(pattern, content)
            matches = matches[0].split()
            matches = [float(match) for match in matches]

            # find etot #

            for line in content.splitlines():
                if re.match(r"\(\s*\d+,\s*\d+\s*\)\s*\|\s*-?\d+\.\d+", line):
                    indices, value = line.split("|")
                    i, j = map(int, re.findall(r"\d+", indices))
                    try:
                        value = float(value.strip())
                        if not (
                            (i == 1 and j == 2)
                            or (i == 2 and j == 1)
                            or (i == 2 and j == 2)
                        ):
                            # res[(i, j)]["etot"] = value
                            res_tmp.append([i, j, value])
                    except:
                        continue

            # print(res_tmp)

            assert len(res_tmp) == len(matches) * raw_data_factor

            res = {}

            for idxqmin, qmin in enumerate(matches):
                res[qmin] = {}

                for i in range(
                    idxqmin * raw_data_factor, idxqmin * raw_data_factor + 5
                ):
                    res[qmin][(res_tmp[i][0], res_tmp[i][1])] = {
                        "etot": 0.0,
                        "ept": res_tmp[i][2],
                        "etot_lcua": 0.0,
                        "etot_full": 0.0,
                    }

                for i in range(
                    idxqmin * raw_data_factor + 5, idxqmin * raw_data_factor + 10
                ):
                    res[qmin][(res_tmp[i][0], res_tmp[i][1])]["etot"] = res_tmp[i][2]

                if with_LCUA:
                    for i in range(
                        idxqmin * raw_data_factor + 10, idxqmin * raw_data_factor + 15
                    ):
                        res[qmin][(res_tmp[i][0], res_tmp[i][1])]["etot_lcua"] = (
                            res_tmp[i][2]
                        )

                if with_full_gFock:
                    for i in range(
                        idxqmin * raw_data_factor + 15, idxqmin * raw_data_factor + 20
                    ):
                        res[qmin][(res_tmp[i][0], res_tmp[i][1])]["etot_full"] = (
                            res_tmp[i][2]
                        )

            return res

        else:
            print(f"Error: File {filename} does not contain the required string.")


def print_out_results(
    nevpt2s: dict, with_LCUA: bool = False, with_full_gFock: bool = False
):
    for qmin, data in nevpt2s.items():
        print(f"Qmin: {qmin}")
        # print ept and etot as a table
        if with_full_gFock:
            table_data = [
                [
                    i,
                    j,
                    value["ept"],
                    value["etot"],
                    value["etot_lcua"],
                    value["etot_full"],
                ]
                for (i, j), value in data.items()
            ]
            headers = ["ncore", "nvirt", "ept", "etot", "etot_lcua", "etot_full"]
        else:
            if with_LCUA:
                table_data = [
                    [i, j, value["ept"], value["etot"], value["etot_lcua"]]
                    for (i, j), value in data.items()
                ]
                headers = ["ncore", "nvirt", "ept", "etot", "etot_lcua"]
            else:
                table_data = [
                    [i, j, value["ept"], value["etot"]]
                    for (i, j), value in data.items()
                ]
                headers = ["ncore", "nvirt", "ept", "etot"]
        print(tabulate(table_data, headers, tablefmt="grid", floatfmt="15.12f"))


############################################
## case II multi qmin
############################################

############################################
## FUll Driver for mr-enpt2/mr-nevpt2d
############################################

from pyscf_util.misc.picture import get_extra_res, draw_extra_pic


def analysis_mrpt2(
    ENPT2_RES,
    MRENPT2_RES,
    NEVPT2_RES,
    TASK,
    CMIN,
    print_pic=False,
    print_table=True,
    NPT=5,
    NLAST_REMOVE_LARGE=None,
    NLAST_REMOVE_SMALL=None,
):

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

    for mole in TASK:
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
        for cmin in CMIN:

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

        # print(DATA_EXTRA)

        DATA_PRINT2 = {}
        DATA_PRINT3 = {}

        data_large_linear = ["large-linear", "energy"]
        data_large_linear_error = ["large-linear", "error"]
        data_large_quadratic = ["large-quadratic", "energy"]
        data_large_quadratic_error = ["large-quadratic", "error"]

        data_small_linear = ["small-linear", "energy"]
        data_small_linear_error = ["small-linear", "error"]
        data_small_quadratic = ["small-quadratic", "energy"]
        data_small_quadratic_error = ["small-quadratic", "error"]

        for key in subspace_order_print:
            DATA_PRINT2[key] = get_extra_res(
                DATA_EXTRA[key]["ept"], DATA_EXTRA[key]["etot"], NPT, NLAST_REMOVE_SMALL
            )
            DATA_PRINT3[key] = get_extra_res(
                DATA_EXTRA[key]["ept"], DATA_EXTRA[key]["etot"], NPT, NLAST_REMOVE_LARGE
            )

            data_large_linear.append(DATA_PRINT3[key]["linear_extra"])
            data_large_linear_error.append(DATA_PRINT3[key]["linear_error"])
            data_large_quadratic.append(DATA_PRINT3[key]["quadratic_extra"])
            data_large_quadratic_error.append(DATA_PRINT3[key]["quadratic_error"])

            data_small_linear.append(DATA_PRINT2[key]["linear_extra"])
            data_small_linear_error.append(DATA_PRINT2[key]["linear_error"])
            data_small_quadratic.append(DATA_PRINT2[key]["quadratic_extra"])
            data_small_quadratic_error.append(DATA_PRINT2[key]["quadratic_error"])

        data_print.append(data_large_linear)
        data_print.append(data_large_quadratic)
        data_print.append(data_small_linear)
        data_print.append(data_small_quadratic)

        data_print.append(data_large_linear_error)
        data_print.append(data_large_quadratic_error)
        data_print.append(data_small_linear_error)
        data_print.append(data_small_quadratic_error)

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

            print("extra with large cmin")
            draw_extra_pic(DATA_PRINT3, 2, 4, subspace_order_print, 24, 9)
            draw_extra_pic(
                DATA_PRINT3, 2, 4, subspace_order_print, 24, 9, use_quadratic=True
            )


def analysis_mrpt2_2(
    ENPT2_RES,
    MRENPT2_RES,
    NEVPT2_RES,
    TASK,
    CMIN,
    print_pic=False,
    print_table=True,
    NPT=5,
    NLAST_REMOVE_LARGE=None,
    NLAST_REMOVE_SMALL=None,
):

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
        for cmin in CMIN:

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

        # print(DATA_EXTRA)

        DATA_PRINT2 = {}
        DATA_PRINT3 = {}

        data_large_linear = ["large-linear", "energy"]
        data_large_linear_error = ["large-linear", "error"]
        data_large_quadratic = ["large-quadratic", "energy"]
        data_large_quadratic_error = ["large-quadratic", "error"]

        data_small_linear = ["small-linear", "energy"]
        data_small_linear_error = ["small-linear", "error"]
        data_small_quadratic = ["small-quadratic", "energy"]
        data_small_quadratic_error = ["small-quadratic", "error"]

        for key in subspace_order_print:
            DATA_PRINT2[key] = get_extra_res(
                DATA_EXTRA[key]["ept"], DATA_EXTRA[key]["etot"], NPT, NLAST_REMOVE_SMALL
            )
            DATA_PRINT3[key] = get_extra_res(
                DATA_EXTRA[key]["ept"], DATA_EXTRA[key]["etot"], NPT, NLAST_REMOVE_LARGE
            )

            data_large_linear.append(DATA_PRINT3[key]["linear_extra"])
            data_large_linear_error.append(DATA_PRINT3[key]["linear_error"])
            data_large_quadratic.append(DATA_PRINT3[key]["quadratic_extra"])
            data_large_quadratic_error.append(DATA_PRINT3[key]["quadratic_error"])

            data_small_linear.append(DATA_PRINT2[key]["linear_extra"])
            data_small_linear_error.append(DATA_PRINT2[key]["linear_error"])
            data_small_quadratic.append(DATA_PRINT2[key]["quadratic_extra"])
            data_small_quadratic_error.append(DATA_PRINT2[key]["quadratic_error"])

        data_print.append(data_large_linear)
        data_print.append(data_large_quadratic)
        data_print.append(data_small_linear)
        data_print.append(data_small_quadratic)

        data_print.append(data_large_linear_error)
        data_print.append(data_large_quadratic_error)
        data_print.append(data_small_linear_error)
        data_print.append(data_small_quadratic_error)

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

            print("extra with large cmin")
            draw_extra_pic(DATA_PRINT3, 2, 2, subspace_order_print, 12, 9)
            draw_extra_pic(
                DATA_PRINT3, 2, 2, subspace_order_print, 12, 9, use_quadratic=True
            )


def extra_nevpt2s(
    DataSet,
    QMIN,
    etot_keyname="etot",
    NPT=5,
    NLAST_REMOVE=None,
    sum=False,
    draw_pic=True,
):

    subspace_order = [(0, 1), (0, 2), (1, 0), (1, 1), (2, 0)]
    subspace_order_2_key = {
        (0, 1): "r",
        (0, 2): "rs",
        (1, 0): "i",
        (1, 1): "ir",
        (2, 0): "ij",
    }

    Res = {}
    for key in subspace_order:
        Ept = [DataSet[x][key]["ept"] for x in QMIN]
        Etot = [DataSet[x][key][etot_keyname] for x in QMIN]

        Res[key] = get_extra_res(Ept, Etot, NPT, NLAST_REMOVE)

    if sum:
        Ept = Res[(0, 1)]["ept"]
        Etot = Res[(0, 1)]["etot"]
        for loc in range(len(QMIN)):
            for key in subspace_order[1:]:
                Ept[loc] += Res[key]["ept"][loc]
                Etot[loc] += Res[key]["etot"][loc]
        Res = get_extra_res(Ept, Etot, NPT, NLAST_REMOVE)

        if draw_pic:
            DATA_DRAW = {}
            DATA_DRAW["ALL"] = Res
            draw_extra_pic(DATA_DRAW, 1, 1, None, 12, 9)
            draw_extra_pic(DATA_DRAW, 1, 1, None, 12, 9, use_quadratic=True)
    else:

        if draw_pic:
            DATA_DRAW = {}
            for key in subspace_order:
                DATA_DRAW[subspace_order_2_key[key]] = Res[key]
            draw_extra_pic(DATA_DRAW, 2, 3, None, 16, 16)
            draw_extra_pic(DATA_DRAW, 2, 3, None, 16, 16, use_quadratic=True)

    return Res


def collect_nevpt2s(
    DataSet,
    TASK,
    CMIN,
    QMIN,
    etot_keyname="etot",
    NPT=5,
    NLAST_REMOVE=None,
    sum=False,
    draw_pic=False,
):

    LinearExtraRes = {}
    LinearExtraErr = {}
    QuadExtraRes = {}
    QuadExtraErr = {}

    subspace_order = [(0, 1), (0, 2), (1, 0), (1, 1), (2, 0)]

    for mole in TASK:
        for cmin in CMIN:
            if draw_pic:
                print(mole, cmin)
            res = extra_nevpt2s(
                DataSet[(mole, cmin)],
                QMIN,
                etot_keyname,
                NPT,
                NLAST_REMOVE,
                sum,
                draw_pic,
            )
            # print(res)
            KEY = (mole, cmin)
            if not sum:
                LinearExtraRes[KEY] = {}
                LinearExtraErr[KEY] = {}
                QuadExtraRes[KEY] = {}
                QuadExtraErr[KEY] = {}
                for key in subspace_order:
                    LinearExtraRes[KEY][key] = res[key]["weighted_linear_extra"]
                    LinearExtraErr[KEY][key] = res[key]["weighted_linear_error"]
                    QuadExtraRes[KEY][key] = res[key]["weighted_quadratic_extra"]
                    QuadExtraErr[KEY][key] = res[key]["weighted_quadratic_error"]
            else:
                LinearExtraRes[KEY] = res["weighted_linear_extra"]
                LinearExtraErr[KEY] = res["weighted_linear_error"]
                QuadExtraRes[KEY] = res["weighted_quadratic_extra"]
                QuadExtraErr[KEY] = res["weighted_quadratic_error"]

    return LinearExtraRes, LinearExtraErr, QuadExtraRes, QuadExtraErr


if __name__ == "__main__":

    filename = "mr_dyall.out"
    matrix_elements = extract_mrenpt2_res(filename)
    print(matrix_elements)

    filename = "mr_enpt2.out"
    matrix_elements = extract_mrenpt2_res(filename)
    print(matrix_elements)

    filename = "mr_sel.out1"
    nevpt2s = extract_nevpt2s_old_type(filename)
    print(nevpt2s)

    filename = "mr_sel.out2"
    nevpt2s = extract_nevpt2s_new_type(filename)
    print(nevpt2s)

    filename = "mr_sel.out3"
    nevpt2s = extract_nevpt2s_new_type(filename, True)
    print_out_results(nevpt2s, True)
    # print(nevpt2s)

    filename = "mr_sel.out4"
    nevpt2s = extract_nevpt2s_new_type(filename, True, True)
    # print(nevpt2s)
    print_out_results(nevpt2s, True, True)
