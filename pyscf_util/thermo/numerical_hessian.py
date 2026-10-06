import numpy as np


# ============================================================
# 工具函数
# ============================================================
def numerical_grad_hess(f, x, step=1e-4, verbose=False):
    """数值梯度与 Hessian（中心差分）。"""
    x = np.asarray(x, dtype=float)
    n = x.size
    h = step

    grad = np.zeros(n)
    for i in range(n):
        if (verbose):
            print(f"grad component {i+1}/{n} done")
        xp = x.copy()
        xm = x.copy()
        xp[i] += h
        xm[i] -= h
        grad[i] = (f(xp) - f(xm)) / (2.0 * h)

    hess = np.zeros((n, n))
    f0 = f(x)
    for i in range(n):
        xp = x.copy()
        xm = x.copy()
        xp[i] += h
        xm[i] -= h
        hess[i, i] = (f(xp) - 2.0 * f0 + f(xm)) / (h * h)

    for i in range(n):
        if verbose:
            print(f"hess row {i+1}/{n} computing")
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


def mass_weighted_frequencies(hess_flat, masses_amu):
    """质量加权 Hessian 对角化 -> 频率 (cm^-1)。"""
    m3 = np.repeat(masses_amu, 3)
    sqrt_m = np.sqrt(m3)
    H_mw = hess_flat / np.outer(sqrt_m, sqrt_m)
    H_mw = 0.5 * (H_mw + H_mw.T)
    eigvals, modes = np.linalg.eigh(H_mw)
    CONV = 5140.4871  # sqrt(Hartree/(Bohr^2·amu)) -> cm^-1
    freqs = np.sign(eigvals) * np.sqrt(np.abs(eigvals)) * CONV
    return freqs, modes


# ============================================================
# 核心：抽象为一个函数
# ============================================================
def compare_numerical_analytic_hess(
    mol,
    mf,
    hess_num,
    grad_num=None,
    step_desc=None,
    freq_zero_thresh=20.0,
    verbose=True,
):
    """
    对比数值 Hessian 与解析 Hessian（含频率）。

    参数
    ----
    mol : pyscf.gto.Mole
        分子对象（坐标必须和 hess_num 一致）。
    mf : pyscf.scf.hf.SCF
        收敛好的 SCF 对象（用于取解析梯度与 Hessian）。
    hess_num : ndarray (3N, 3N)
        数值 Hessian，单位 Hartree/Bohr^2。
        坐标顺序 [atom0_x, atom0_y, atom0_z, atom1_x, ...]。
    grad_num : ndarray (3N,) or (natm,3) or None
        数值梯度，单位 Hartree/Bohr。None 时跳过梯度对比。
    step_desc : str or None
        步长描述，仅用于打印。
    freq_zero_thresh : float
        低于该值的频率视为平动/转动模，单位 cm^-1。默认 20。
    verbose : bool
        是否打印详细信息。

    返回
    ----
    result : dict
        包含所有中间量与对比结果：
        'grad_ana', 'grad_num', 'hess_ana', 'hess_num',
        'freqs_ana', 'freqs_num',
        'freqs_ana_real', 'freqs_num_real',
        'freq_diff_max', 'freq_diff_rms',
        'hess_diff_max', 'grad_diff_max'
    """
    natm = mol.natm
    ndim = 3 * natm
    masses = np.array(mol.atom_mass_list(isotope_avg=True))

    # ---------- 解析梯度 ----------
    grad_ana = mf.nuc_grad_method().kernel()  # (natm,3)
    grad_ana_flat = grad_ana.reshape(-1)

    # ---------- 解析 Hessian ----------
    hess_ana = mf.Hessian().kernel()  # (natm,natm,3,3)
    hess_ana_flat = hess_ana.transpose(0, 2, 1, 3).reshape(ndim, ndim)

    # ---------- 数值梯度处理 ----------
    if grad_num is not None:
        grad_num_flat = np.asarray(grad_num).reshape(-1)
    else:
        grad_num_flat = None

    # ============================================================
    # 打印工具
    # ============================================================
    def _mat_diff(name, A, B):
        diff = A - B
        rms = np.sqrt(np.mean(diff**2))
        maxd = np.max(np.abs(diff))
        maxref = np.max(np.abs(B))
        rel = maxd / max(maxref, 1e-30)
        print(f"\n=== {name} ===")
        print(f"  shape        : {A.shape}")
        print(f"  max |diff|   : {maxd:.3e}")
        print(f"  rms  diff    : {rms:.3e}")
        print(f"  max |ref|    : {maxref:.3e}")
        print(f"  max rel diff : {rel:.3e}")
        return maxd, rms

    def _print_freqs(freqs, label):
        print(f"\n=== {label} ===")
        for i, f in enumerate(freqs):
            tag = ""
            if abs(f) < freq_zero_thresh:
                tag = "  <- 平动/转动"
            elif f < 0:
                tag = "  <- 虚频"
            print(f"    {i:3d}  {f:14.4f}{tag}")

    # ============================================================
    # 梯度 / Hessian 对比
    # ============================================================
    if step_desc is not None:
        print(f"\n[step = {step_desc}]")

    grad_diff_max = None
    if grad_num_flat is not None:
        grad_diff_max, _ = _mat_diff(
            "Gradient: numerical vs analytic", grad_num_flat, grad_ana_flat
        )
        if verbose:
            print("\n解析梯度 (Hartree/Bohr):")
            print(grad_ana)
            print("\n数值梯度 (Hartree/Bohr):")
            print(grad_num_flat.reshape(natm, 3))

    hess_diff_max, _ = _mat_diff(
        "Hessian: numerical vs analytic", hess_num, hess_ana_flat
    )
    if verbose:
        np.set_printoptions(precision=8, suppress=True)
        print("\n解析 Hessian (Hartree/Bohr^2):")
        print(hess_ana_flat)
        print("\n数值 Hessian (Hartree/Bohr^2):")
        print(hess_num)

    # ============================================================
    # 频率对比
    # ============================================================
    freqs_ana, _ = mass_weighted_frequencies(hess_ana_flat, masses)
    freqs_num, _ = mass_weighted_frequencies(hess_num, masses)

    idx_a = np.argsort(freqs_ana)
    idx_n = np.argsort(freqs_num)
    freqs_ana_sorted = freqs_ana[idx_a]
    freqs_num_sorted = freqs_num[idx_n]

    if verbose:
        _print_freqs(freqs_ana_sorted, "Analytic Hessian frequencies")
        _print_freqs(freqs_num_sorted, "Numerical Hessian frequencies")

    mask_a = np.abs(freqs_ana_sorted) > freq_zero_thresh
    mask_n = np.abs(freqs_num_sorted) > freq_zero_thresh
    fa = freqs_ana_sorted[mask_a]
    fn = freqs_num_sorted[mask_n]

    freq_diff_max = None
    freq_diff_rms = None

    print("\n=== Frequency comparison (cm^-1) ===")
    print(f"  解析非零模数: {fa.size}, 数值非零模数: {fn.size}")
    if fa.size == fn.size and fa.size > 0:
        diff = fa - fn
        print(f"  {'#':>3} {'analytic':>12} {'numerical':>12} {'diff':>12}")
        for i, (a, n, d) in enumerate(zip(fa, fn, diff)):
            print(f"  {i:>3} {a:>12.4f} {n:>12.4f} {d:>12.4f}")
        freq_diff_max = float(np.max(np.abs(diff)))
        freq_diff_rms = float(np.sqrt(np.mean(diff**2)))
        print(f"\n  max |Δfreq| = {freq_diff_max:.4f} cm^-1")
        print(f"  rms  Δfreq  = {freq_diff_rms:.4f} cm^-1")
    else:
        print("  非零模数不一致，跳过逐项对比。")

    print("\n平动/转动模式的近零频率 (cm^-1):")
    print(
        "  analytic :", freqs_ana_sorted[np.abs(freqs_ana_sorted) <= freq_zero_thresh]
    )
    print(
        "  numerical:", freqs_num_sorted[np.abs(freqs_num_sorted) <= freq_zero_thresh]
    )

    return {
        "grad_ana": grad_ana,
        "grad_num": grad_num_flat,
        "hess_ana": hess_ana_flat,
        "hess_num": hess_num,
        "freqs_ana": freqs_ana_sorted,
        "freqs_num": freqs_num_sorted,
        "freqs_ana_real": fa,
        "freqs_num_real": fn,
        "grad_diff_max": grad_diff_max,
        "hess_diff_max": hess_diff_max,
        "freq_diff_max": freq_diff_max,
        "freq_diff_rms": freq_diff_rms,
    }


# ============================================================
# PySCF 使用示例
# ============================================================
if __name__ == "__main__":
    from pyscf import gto, scf

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
    x0 = mol.atom_coords(unit="Bohr").ravel()

    # ---------- 解析梯度 / Hessian ----------
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-10
    mf.max_cycle = 200
    mf.kernel()
    assert mf.converged

    # ---------- 数值梯度 / Hessian ----------
    def energy(x):
        mol.set_geom_(x.reshape(-1, 3), unit="Bohr")
        m = scf.RHF(mol)
        m.conv_tol = 1e-12
        m.conv_tol_grad = 1e-10
        m.max_cycle = 200
        e = m.kernel()
        if not m.converged:
            print("WARNING: SCF not converged during numerical diff!")
        return e

    step = 1e-3
    grad_num, hess_num = numerical_grad_hess(energy, x0, step=step, verbose=True)

    # ⚠️ 数值差分后 mol 的几何被改成了最后一次扰动的坐标，
    #    需要恢复初始几何，再算解析量。
    mol.set_geom_(x0.reshape(-1, 3), unit="Bohr")

    # ---------- 一步调用对比 ----------

    result = compare_numerical_analytic_hess(
        mol=mol,
        mf=mf,
        hess_num=hess_num,
        grad_num=grad_num,
        step_desc=f"{step:.1e} Bohr",
        freq_zero_thresh=20.0,
        verbose=False,
    )

    # 需要结构化结果时:
    # print("max Δfreq =", result["freq_diff_max"])
