import numpy as np
from pyscf import gto, scf


def numerical_grad_hess(f, x, step=1e-4, verbose=False):
    """
    对标量函数 f(x) 计算数值梯度和 Hessian。

    参数
    ----
    f : callable
        输入一维数组 x，返回标量 float。
    x : array_like
        一维数组，长度 n。
    step : float
        有限差分步长。能量对核坐标通常取 1e-4 ~ 1e-5 Bohr。
    verbose : bool
        是否打印进度。

    返回
    ----
    grad : ndarray, shape (n,)
        数值梯度，中心差分。
    hess : ndarray, shape (n, n)
        数值 Hessian，中心差分，对称矩阵。
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    h = step

    if verbose:
        print(f"n = {n}, step = {h}")

    # ---------- 梯度 ----------
    grad = np.zeros(n)
    for i in range(n):
        xp = x.copy()
        xm = x.copy()
        xp[i] += h
        xm[i] -= h
        grad[i] = (f(xp) - f(xm)) / (2.0 * h)
        if verbose:
            print(f"grad {i+1}/{n} done")

    # ---------- Hessian ----------
    hess = np.zeros((n, n))
    f0 = f(x)

    # 对角线：二阶中心差分
    for i in range(n):
        xp = x.copy()
        xm = x.copy()
        xp[i] += h
        xm[i] -= h
        hess[i, i] = (f(xp) - 2.0 * f0 + f(xm)) / (h * h)
        if verbose:
            print(f"hess diag {i+1}/{n} done")

    # 非对角线：混合中心差分
    for i in range(n):
        for j in range(i + 1, n):
            xpp = x.copy()
            xpm = x.copy()
            xmp = x.copy()
            xmm = x.copy()

            xpp[i] += h
            xpp[j] += h
            xpm[i] += h
            xpm[j] -= h
            xmp[i] -= h
            xmp[j] += h
            xmm[i] -= h
            xmm[j] -= h

            val = (f(xpp) - f(xpm) - f(xmp) + f(xmm)) / (4.0 * h * h)
            hess[i, j] = val
            hess[j, i] = val
        if verbose:
            print(f"hess row {i+1}/{n} done")

    return grad, hess


def print_matrix_diff(name, A, B):
    """打印两个矩阵的对比信息。"""
    diff = A - B
    print(f"\n=== {name} ===")
    print(f"  shape        : {A.shape}")
    print(f"  max |diff|   : {np.max(np.abs(diff)):.3e}")
    print(f"  rms  diff    : {np.sqrt(np.mean(diff**2)):.3e}")
    print(f"  max |ref|    : {np.max(np.abs(B)):.3e}")
    rel = np.max(np.abs(diff)) / max(np.max(np.abs(B)), 1e-30)
    print(f"  max rel diff : {rel:.3e}")


def mass_weighted_frequencies(hess_flat, masses_amu):
    """
    对 Hessian 做质量加权对角化，返回频率（cm^-1）。

    参数
    ----
    hess_flat : ndarray (3N, 3N)
        Hessian，单位 Hartree/Bohr^2。坐标顺序是
        [atom0_x, atom0_y, atom0_z, atom1_x, ...]。
    masses_amu : ndarray (N,)
        原子质量，单位 amu。

    返回
    ----
    freqs  : ndarray (3N,)   频率，单位 cm^-1（含虚频为负值）
    modes  : ndarray (3N,3N) 简正模式（每一列是一个模式）
    """
    # 每个坐标对应的质量
    m3 = np.repeat(masses_amu, 3)
    # 质量加权：H_mw = M^{-1/2} H M^{-1/2}
    sqrt_m = np.sqrt(m3)
    H_mw = hess_flat / np.outer(sqrt_m, sqrt_m)

    # 对称化（数值误差可能引入很小反对称部分）
    H_mw = 0.5 * (H_mw + H_mw.T)

    eigvals, modes = np.linalg.eigh(H_mw)

    # 频率转换因子：sqrt(Hartree/(Bohr^2*amu)) -> cm^-1
    # sqrt(eig) * 219474.63 (Hartree->cm^-1) / (2*pi*c) 合并成常数
    # 常数为 5140.4871 cm^-1 / sqrt(Hartree/(Bohr^2*amu))
    CONV = 5140.4871
    freqs = np.sign(eigvals) * np.sqrt(np.abs(eigvals)) * CONV
    return freqs, modes


def print_frequencies(freqs, label, n_trans_rot=6):
    """打印频率，标注平动转动（近零）和虚频。"""
    print(f"\n=== {label} ===")
    print("  模式频率 (cm^-1):")
    for i, f in enumerate(freqs):
        tag = ""
        if abs(f) < 20:
            tag = "  <- 平动/转动"
        elif f < 0:
            tag = "  <- 虚频"
        print(f"    {i:3d}  {f:14.4f}{tag}")


# ============================================================
# PySCF 示例：计算 RHF 的数值梯度与 Hessian
# ============================================================
if __name__ == "__main__":
    # 1. 定义分子
    mol = gto.M(
        atom="""
        O  0.000000  0.000000  0.000000
        H  0.000000  0.000000  0.9572
        H  0.926627  0.000000 -0.239987
        """,
        basis="6-31g",
        unit="Angstrom",
        verbose=0,
    )

    # 初始核坐标，展平为一维向量，单位 Bohr
    x0 = mol.atom_coords(unit="Bohr").ravel()
    natm = mol.natm

    # 原子质量 (amu)。也可以用 mol.atom_mass_list(isotope_avg=True)
    masses = np.array(mol.atom_mass_list(isotope_avg=True))

    # 2. 定义能量函数：输入核坐标向量 x，返回 RHF 能量
    def energy(x):
        # 更新分子几何
        mol.set_geom_(x.reshape(-1, 3), unit="Bohr")
        # 每次都重新做 RHF，保证从当前几何出发
        mf = scf.RHF(mol)
        mf.conv_tol = 1e-14
        mf.conv_tol_grad = 1e-10
        mf.max_cycle = 200
        e = mf.kernel()
        return e

    # ---------- 解析梯度 ----------
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-10
    mf.max_cycle = 200
    mf.kernel()

    # 解析梯度：单位 Hartree/Bohr，shape (natm, 3)
    grad_ana = mf.nuc_grad_method().kernel()
    grad_ana_flat = grad_ana.reshape(-1)

    # 解析 Hessian：单位 Hartree/Bohr^2，shape (natm, natm, 3, 3)
    hess_ana = mf.Hessian().kernel()
    # 展平为 (3N, 3N)，顺序与 atom_coords.ravel() 对应
    hess_ana_flat = hess_ana.transpose(0, 2, 1, 3).reshape(3 * natm, 3 * natm)

    # 3. 计算数值梯度和 Hessian
    grad, hess = numerical_grad_hess(energy, x0, step=1e-3, verbose=True)

    np.set_printoptions(precision=8, suppress=True)

    print("\n=== Numerical Gradient (Hartree/Bohr) ===")
    print(grad.reshape(natm, 3))

    print("\n=== Analytic Gradient (Hartree/Bohr) ===")
    print(grad_ana.reshape(natm, 3))

    print("\n=== Numerical Hessian (Hartree/Bohr^2) ===")
    print(hess)

    print("\n=== Analytic Hessian (Hartree/Bohr^2) ===")
    print(hess_ana_flat)

    # 可选：重塑为 PySCF 常见形状 (natm, natm, 3, 3)
    hess_reshaped = hess.reshape(natm, 3, natm, 3).transpose(0, 2, 1, 3)
    print("\nHessian shape (natm, natm, 3, 3):", hess_reshaped.shape)
    # print(hess_reshaped)

    print_matrix_diff(
        "Gradient: numerical vs analytic", grad.reshape(-1), grad_ana_flat
    )
    print_matrix_diff("Hessian: numerical vs analytic", hess, hess_ana_flat)

    # ---------- 对角化 -> 频率 ----------
    freqs_ana, modes_ana = mass_weighted_frequencies(hess_ana_flat, masses)
    freqs_num, modes_num = mass_weighted_frequencies(hess, masses)

    # 排序（eigh 返回升序，虚频是负数，会排在最前面）
    idx_a = np.argsort(freqs_ana)
    idx_n = np.argsort(freqs_num)

    print_frequencies(freqs_ana[idx_a], "Analytic Hessian frequencies")
    print_frequencies(freqs_num[idx_n], "Numerical Hessian frequencies")

    # ---------- 频率对比 ----------
    # 只比较非平动/转动模式 (|freq| > 20 cm^-1)
    mask_a = np.abs(freqs_ana) > 20
    mask_n = np.abs(freqs_num) > 20
    fa = freqs_ana[mask_a]
    fn = freqs_num[mask_n]

    print("\n=== Frequency comparison (cm^-1) ===")
    print(f"  解析非零模数: {fa.size}, 数值非零模数: {fn.size}")
    if fa.size == fn.size:
        diff = fa - fn
        print(f"  {'#':>3} {'analytic':>12} {'numerical':>12} {'diff':>12}")
        for i, (a, n, d) in enumerate(zip(fa, fn, diff)):
            print(f"  {i:>3} {a:>12.4f} {n:>12.4f} {d:>12.4f}")
        print(f"\n  max |Δfreq| = {np.max(np.abs(diff)):.4f} cm^-1")
        print(f"  rms  Δfreq  = {np.sqrt(np.mean(diff**2)):.4f} cm^-1")
    else:
        print("  非零模数不一致，跳过逐项对比。")

    # 可选：把平动/转动模数打印出来
    print("\n平动/转动模式的近零特征值 (cm^-1):")
    print("  analytic :", freqs_ana[np.abs(freqs_ana) <= 20])
    print("  numerical:", freqs_num[np.abs(freqs_num) <= 20])
