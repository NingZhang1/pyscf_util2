import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


class BasisExtrapolation:
    """
    能量外推类，包含HF能量的指数外推和关联能的X^-3外推
    """

    def __init__(self, basis_sets, total_energies, hf_energies):
        """
        初始化

        参数:
        basis_sets: 基组列表，如['DZ', 'TZ', 'QZ', '5Z']
        total_energies: 总能量列表
        hf_energies: HF能量列表
        """
        self.basis_sets = basis_sets
        self.total_energies = np.array(total_energies)
        self.hf_energies = np.array(hf_energies)

        # 检查输入数据
        if len(basis_sets) != len(total_energies) or len(basis_sets) != len(
            hf_energies
        ):
            raise ValueError("所有输入列表的长度必须相同")

        # 基组对应的数值表示
        self.basis_values = {"DZ": 2, "TZ": 3, "QZ": 4, "5Z": 5}
        self.X = np.array(
            [self.basis_values.get(basis, i + 2) for i, basis in enumerate(basis_sets)]
        )

    def calculate_correlation_energies(self):
        """计算关联能"""
        return self.total_energies - self.hf_energies

    def exponential_func(self, x, a, b, c):
        """指数函数用于HF能量外推: E = a + b * exp(-c*x)"""
        return a + b * np.exp(-c * x)

    def extrapolate_hf_energy(self):
        """
        对HF能量进行指数外推

        返回:
        hf_infinity: 外推到无限基组的HF能量
        popt: 拟合参数
        """
        # 设置初始猜测参数
        p0 = [
            min(self.hf_energies),
            max(self.hf_energies) - min(self.hf_energies),
            np.random.rand() * 0.5 + 0.5,
        ]

        try:
            # 使用curve_fit进行非线性最小二乘拟合
            popt, pcov = curve_fit(
                self.exponential_func,
                self.X[-3:],
                self.hf_energies[-3:],
                p0=p0,
                maxfev=5000,
            )
            hf_infinity = popt[0]  # a参数对应无限基组的HF能量
            # print(popt)
            # print(pcov)

            # 生成外推曲线
            x_fit = np.linspace(min(self.X), max(self.X) + 2, 100)
            y_fit = self.exponential_func(x_fit, *popt)

            return hf_infinity, popt, x_fit, y_fit

        except Exception as e:
            print(f"指数拟合失败: {e}")
            # 如果非线性拟合失败，尝试使用线性化的指数拟合
            return self.extrapolate_hf_energy_linearized()

    def extrapolate_hf_energy_linearized(self):
        """
        线性化的指数拟合方法
        """
        # 寻找HF能量的极限估计
        hf_diff = np.diff(self.hf_energies)
        hf_infinity_est = (
            self.hf_energies[-1] + hf_diff[-1] ** 2 / (hf_diff[-1] - hf_diff[-2])
            if len(hf_diff) > 1
            else self.hf_energies[-1]
        )

        # 简单指数衰减拟合
        popt = [hf_infinity_est, self.hf_energies[0] - hf_infinity_est, 0.5]
        x_fit = np.linspace(min(self.X), max(self.X) + 2, 100)
        y_fit = self.exponential_func(x_fit, *popt)

        return hf_infinity_est, popt, x_fit, y_fit

    def extrapolate_correlation_energy(self, correlation_energies):
        """
        对关联能进行X^-3外推（使用最后两个点：QZ和5Z）

        参数:
        correlation_energies: 关联能数组

        返回:
        corr_infinity: 外推到无限基组的关联能
        C: 拟合参数
        """
        # 找到QZ和5Z对应的索引
        qz_idx = None
        z5_idx = None

        for i, basis in enumerate(self.basis_sets):
            if basis == "QZ":
                qz_idx = i
            elif basis == "5Z":
                z5_idx = i

        if qz_idx is None or z5_idx is None:
            # 如果没有明确的QZ和5Z标签，使用最后两个点
            print("警告: 未找到QZ和5Z标签，使用最后两个数据点进行外推")
            qz_idx = -2
            z5_idx = -1

        # 获取最后两个点的数据
        X_last_two = self.X[[qz_idx, z5_idx]]
        E_corr_last_two = correlation_energies[[qz_idx, z5_idx]]

        # 使用X^-3外推公式: E_corr = E_corr_inf + C * X^{-3}
        # 解方程:
        # E1 = E_inf + C * X1^{-3}
        # E2 = E_inf + C * X2^{-3}

        X1_inv3 = float(X_last_two[0]) ** (-3)
        X2_inv3 = float(X_last_two[1]) ** (-3)
        E1 = E_corr_last_two[0]
        E2 = E_corr_last_two[1]

        # 解线性方程组
        # C = (E1 - E2) / (X1^{-3} - X2^{-3})
        # E_inf = E1 - C * X1^{-3}

        C = (E1 - E2) / (X1_inv3 - X2_inv3)
        corr_infinity = E1 - C * X1_inv3

        return corr_infinity, C, qz_idx, z5_idx

    def perform_extrapolation(self):
        """
        执行完整的外推过程

        返回:
        results: 包含所有结果的字典
        """
        # 计算关联能
        correlation_energies = self.calculate_correlation_energies()

        # HF能量指数外推
        hf_infinity, hf_params, hf_x_fit, hf_y_fit = self.extrapolate_hf_energy()

        # 关联能X^-3外推
        corr_infinity, C, qz_idx, z5_idx = self.extrapolate_correlation_energy(
            correlation_energies
        )

        # 总能量外推结果
        total_infinity = hf_infinity + corr_infinity

        # 整理结果
        results = {
            "hf_infinity": hf_infinity,
            "hf_params": hf_params,
            "corr_infinity": corr_infinity,
            "C_parameter": C,
            "total_infinity": total_infinity,
            "correlation_energies": correlation_energies,
            "hf_x_fit": hf_x_fit,
            "hf_y_fit": hf_y_fit,
            "qz_idx": qz_idx,
            "z5_idx": z5_idx,
        }

        return results

    def plot_results(self, results):
        """
        绘制外推结果图
        """
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        # 1. HF能量外推图
        ax1 = axes[0, 0]
        ax1.plot(self.X, self.hf_energies, "bo-", label="HF Energy", markersize=8)
        ax1.plot(
            results["hf_x_fit"], results["hf_y_fit"], "r--", label="Exponential Fit"
        )
        ax1.axhline(
            y=results["hf_infinity"],
            color="g",
            linestyle=":",
            label=f'HF∞ = {results["hf_infinity"]:.6f}',
        )
        ax1.set_xlabel("Basis Set (X)")
        ax1.set_ylabel("HF Energy")
        ax1.set_title("HF Energy Exponential Extrapolation")
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # 标记基组点
        for i, (x, basis) in enumerate(zip(self.X, self.basis_sets)):
            ax1.annotate(
                basis,
                (x, self.hf_energies[i]),
                textcoords="offset points",
                xytext=(0, 10),
                ha="center",
            )

        # 2. 关联能外推图
        ax2 = axes[0, 1]
        correlation_energies = results["correlation_energies"]

        # 绘制所有关联能点
        ax2.plot(
            self.X,
            correlation_energies,
            "bo-",
            label="Correlation Energy",
            markersize=8,
        )

        # 特别标记最后两个点
        qz_idx, z5_idx = results["qz_idx"], results["z5_idx"]
        ax2.plot(
            self.X[qz_idx],
            correlation_energies[qz_idx],
            "ro",
            markersize=10,
            label="QZ (for extrapolation)",
        )
        ax2.plot(
            self.X[z5_idx],
            correlation_energies[z5_idx],
            "ro",
            markersize=10,
            label="5Z (for extrapolation)",
        )

        # 绘制外推线
        x_extrap = np.linspace(min(self.X), max(self.X) + 3, 100)
        y_extrap = results["corr_infinity"] + results["C_parameter"] * (
            x_extrap ** (-3)
        )
        ax2.plot(x_extrap, y_extrap, "r--", label="X^-3 Extrapolation")

        ax2.axhline(
            y=results["corr_infinity"],
            color="g",
            linestyle=":",
            label=f'Corr∞ = {results["corr_infinity"]:.6f}',
        )
        ax2.set_xlabel("Basis Set (X)")
        ax2.set_ylabel("Correlation Energy")
        ax2.set_title("Correlation Energy X^-3 Extrapolation")
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        # 3. 总能量图
        ax3 = axes[1, 0]
        ax3.plot(self.X, self.total_energies, "bo-", label="Total Energy", markersize=8)
        ax3.axhline(
            y=results["total_infinity"],
            color="r",
            linestyle="--",
            label=f'Total∞ = {results["total_infinity"]:.6f}',
        )
        ax3.set_xlabel("Basis Set (X)")
        ax3.set_ylabel("Total Energy")
        ax3.set_title("Total Energy with Extrapolated Limit")
        ax3.legend()
        ax3.grid(True, alpha=0.3)

        # 4. 残差图
        ax4 = axes[1, 1]
        hf_fit_values = results["hf_params"][0] + results["hf_params"][1] * np.exp(
            -results["hf_params"][2] * self.X
        )
        hf_residuals = self.hf_energies - hf_fit_values

        ax4.plot(self.X, hf_residuals, "bo-", markersize=8)
        ax4.axhline(y=0, color="r", linestyle="--")
        ax4.set_xlabel("Basis Set (X)")
        ax4.set_ylabel("Residuals")
        ax4.set_title("HF Energy Fit Residuals")
        ax4.grid(True, alpha=0.3)

        plt.tight_layout()
        return fig
