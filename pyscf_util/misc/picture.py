from urllib import robotparser

# import seaborn
# import pandas
import matplotlib.pyplot as plt

# import numpy


def draw_extra_pic(
    x: list,
    y: list,
    legend: list,
    line_prop: list,
    xlabel: str = "$E_{pt}^{(2)}/E_H$",
    ylabel: str = "$E_{tot}/E_H$",
    title="",
    width=16,
    height=9,
    fontsize_xylabel=18,
    fontsize_xytick=18,
    fontsize_title=18,
    fontsize_legend=18,
    save_name=None,
):
    plt.figure(figsize=(width, height))
    for id, x in enumerate(x):
        plt.plot(
            x,
            y[id],
            marker=line_prop[id]["marker"],
            markersize=line_prop[id]["markersize"],
            linewidth=line_prop[id]["linewidth"],
            label=legend[id],
        )
    plt.xlabel(xlabel, fontsize=fontsize_xylabel)
    plt.ylabel(ylabel, fontsize=fontsize_xylabel)
    plt.xticks(fontsize=fontsize_xytick)
    plt.yticks(fontsize=fontsize_xytick)
    plt.title(title, fontsize=fontsize_title)
    plt.legend(fontsize=fontsize_legend)
    if save_name is not None:
        plt.savefig(save_name)
    plt.show()


######################################
# extra res #
######################################

from pyscf_util.Analyzer.iCIPT2.analyzer import *


def get_extra_res(ept, etot, NEXTRA=5, NLAST_REMOVE=None):  #

    if NLAST_REMOVE == 0:
        NLAST_REMOVE = None

    try:

        ### do linear analysis ###

        if NLAST_REMOVE is None:
            # linear #
            _, linear_extra, _, linear_error = LinearRegression_EstimateError(
                ept[-NEXTRA:], etot[-NEXTRA:]
            )
            # weighted linear #
            weighted_linear_extra, weighted_linear_error, slope, intercept = (
                weighted_linear_fit_estimate_error(
                    ept[-NEXTRA:], etot[-NEXTRA:], False, True
                )
            )
            # quad #
            quadratic_extra, quadratic_error = quadratic_fit_estimate_error(
                ept[-NEXTRA:], etot[-NEXTRA:], False, False
            )
            # weighted quad #
            weighted_quadratic_extra, weighted_quadratic_error, a, b, c = (
                weighted_quadratic_fit_estimate_error(
                    ept[-NEXTRA:], etot[-NEXTRA:], False, True
                )
            )

        else:
            # linear #
            _, linear_extra, _, linear_error = LinearRegression_EstimateError(
                ept[-NLAST_REMOVE - NEXTRA : -NLAST_REMOVE],
                etot[-NLAST_REMOVE - NEXTRA : -NLAST_REMOVE],
            )
            # weighted linear #
            weighted_linear_extra, weighted_linear_error, slope, intercept = (
                weighted_linear_fit_estimate_error(
                    ept[-NLAST_REMOVE - NEXTRA : -NLAST_REMOVE],
                    etot[-NLAST_REMOVE - NEXTRA : -NLAST_REMOVE],
                    False,
                    True,
                )
            )
            # quad #
            quadratic_extra, quadratic_error = quadratic_fit_estimate_error(
                ept[-NLAST_REMOVE - NEXTRA : -NLAST_REMOVE],
                etot[-NLAST_REMOVE - NEXTRA : -NLAST_REMOVE],
                False,
                False,
            )
            # weighted quad #
            weighted_quadratic_extra, weighted_quadratic_error, a, b, c = (
                weighted_quadratic_fit_estimate_error(
                    ept[-NLAST_REMOVE - NEXTRA : -NLAST_REMOVE],
                    etot[-NLAST_REMOVE - NEXTRA : -NLAST_REMOVE],
                    False,
                    True,
                )
            )

        ### return res ###

        return {
            "ept": ept,
            "etot": etot,
            "linear_extra": linear_extra,
            "linear_error": linear_error,
            "weighted_linear_extra": weighted_linear_extra,
            "weighted_linear_error": weighted_linear_error,
            "slope": slope,
            "intercept": intercept,
            "a": a,
            "b": b,
            "c": c,
            "quadratic_extra": quadratic_extra,
            "quadratic_error": quadratic_error,
            "weighted_quadratic_extra": weighted_quadratic_extra,
            "weighted_quadratic_error": weighted_quadratic_error,
        }

    except Exception as e:
        print(e)
        return {
            "ept": ept,
            "etot": etot,
            "linear_extra": 0.0,
            "linear_error": 0.0,
            "weighted_linear_extra": 0.0,
            "weighted_linear_error": 0.0,
            "slope": 0.0,
            "intercept": 0.0,
            "a": 0.0,
            "b": 0.0,
            "c": 0.0,
            "quadratic_extra": 0.0,
            "quadratic_error": 0.0,
            "weighted_quadratic_extra": 0.0,
            "weighted_quadratic_error": 0.0,
        }


######################################
# draw picture #
######################################


def draw_extra_pic(
    DATA,
    nfig_x,
    nfig_y,
    subtasks=None,
    figsize_x=16,
    figsize_y=12,
    fig_title=f"iCIPT2 Regression",
    origin_data_label=f"Original Data",
    linear_fit_label=f"Weighted Linear Fit",
    quadratic_fit_label=f"Weighted Quadratic Fit",
    add_err_bar=True,
    use_quadratic=False,
    shown=True,
    save_fig=False,
    fig_path=None,
):
    if subtasks is None:
        subtasks = DATA.keys()

    colors = ["blue", "red", "green", "orange", "purple", "brown"]

    if nfig_x * nfig_y < len(subtasks):
        exit(1)

    ################################################################
    # construct #
    ################################################################

    fig, axes = plt.subplots(nfig_x, nfig_y, figsize=(figsize_x, figsize_y))
    fig.suptitle(
        fig_title,
        fontsize=16,
        fontweight="bold",
    )

    try:
        axes = axes.flatten()
    except Exception as e:
        pass

    ################################################################
    # loop #
    ################################################################

    for idx, taskname in enumerate(subtasks):

        try:
            if idx >= len(axes):
                break
            ax = axes[idx]
        except Exception as e:
            ax = axes
            pass

        data = DATA[taskname]

        ept = data["ept"]
        etot = data["etot"]

        # 绘制原始数据点

        ax.scatter(
            ept,
            etot,
            color=colors[idx % len(colors)],
            s=60,
            alpha=0.7,
            label=origin_data_label,
            zorder=5,
        )

        # 绘制拟合曲线

        if not use_quadratic:

            # 创建拟合直线

            slope = data["slope"]
            intercept = data["intercept"]
            weighted_extra = data["weighted_linear_extra"]
            weighted_error = data["weighted_linear_error"]

            x_fit = np.linspace(min(ept), 0.0, 1000)
            y_fit = slope * x_fit + intercept

            ax.plot(
                x_fit,
                y_fit,
                "--",
                color=colors[idx % len(colors)],
                linewidth=2,
                label=linear_fit_label,
                zorder=4,
            )

        else:

            a = data["a"]
            b = data["b"]
            c = data["c"]
            weighted_extra = data["weighted_quadratic_extra"]
            weighted_error = data["weighted_quadratic_error"]
            x_fit = np.linspace(min(ept), 0.0, 1000)
            y_fit = a * x_fit**2 + b * x_fit + c

            ax.plot(
                x_fit,
                y_fit,
                "--",
                color=colors[idx % len(colors)],
                linewidth=2,
                label=quadratic_fit_label,
                zorder=4,
            )

        # 标记外推点（x=0）

        ax.axvline(x=0, color="gray", linestyle=":", alpha=0.7)
        ax.scatter(
            0,
            weighted_extra,
            color="black",
            s=100,
            marker="*",
            label=f"Extrapolated: {weighted_extra:.6f} ± {weighted_error:.6f}",
            zorder=6,
        )

        # 添加误差条

        if add_err_bar:
            ax.errorbar(
                0,
                weighted_extra,
                yerr=weighted_error,
                fmt="none",
                ecolor="red",
                elinewidth=2,
                capsize=5,
                capthick=2,
            )

        # 设置标签和标题

        ax.set_xlabel(r"$E_{\text{c}}^{(2)}$ / $E_\text{h}$", fontsize=12)
        ax.set_ylabel(r"$E_{\text{tot}}$ / $E_\text{h}$", fontsize=12)
        ax.set_title(f"{taskname}", fontsize=14, fontweight="bold")

        # 添加图例

        ax.legend(loc="best", fontsize=10)

        # 添加网格

        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.subplots_adjust(top=0.93)

    if shown:
        plt.show()

    if save_fig:
        plt.savefig(fig_path)
