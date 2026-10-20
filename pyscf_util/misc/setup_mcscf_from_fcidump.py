import numpy as np
from pyscf import mcscf, tools


def setup_mcscf_from_fcidump(
    fcidump_file,
    ncas,
    nelecas,
    active_orbs=None,
    mc_type="casscf",
    molpro_orbsym=None,
    verbose=0,
):
    """
    从 FCIDUMP 文件读取积分（假定轨道正交），构建 PySCF MCSCF 对象。

    本函数不进行 SCF 迭代：通过设置 ``mf.max_cycle = 0`` 并调用
    ``mf.kernel()``，仅用于构造初始轨道/密度等必要属性。

    Parameters
    ----------
    fcidump_file : str
        FCIDUMP 文件路径。
    ncas : int
        活性空间轨道数。
    nelecas : int or tuple of int
        活性空间的电子数（可以是总电子数，或 (nalpha, nbeta) 元组）。
    active_orbs : list of int, optional
        0-based 轨道索引列表，用于重排分子轨道，使"冻结芯 + 活性"轨道
        排在前 ncore+ncas 位。应给出 **core + active** 轨道的索引
        （顺序不限，函数内部会按你给出的顺序排列）。
        若为 None，则采用 PySCF 默认（按能量排序取最高的 ncas 个为活性）。
    mc_type : {'casscf', 'casci'}
        构建 CASSCF 还是 CASCI 对象。
    molpro_orbsym : str or list, optional
        传递给 ``tools.fcidump.to_scf`` 的 ``molpro_orbsym`` 参数，
        用于处理 Molpro 风格的 ORBSYM 约定。
    verbose : int
        PySCF 输出级别。

    Returns
    -------
    mol : pyscf.gto.Mole
    mf : pyscf.scf.hf.SCF
        已完成 ``kernel()`` 调用（未做 SCF 迭代）的 SCF 对象。
    mc : pyscf.mcscf.CASSCF or pyscf.mcscf.CASCI
    """
    # ---- 1. 读取 FCIDUMP，构建正交轨道模型 ----
    mf = tools.fcidump.to_scf(fcidump_file, molpro_orbsym=molpro_orbsym)
    mf.verbose = verbose
    mol = mf.mol

    # ---- 2. 不迭代 SCF：只调用一次 kernel 构造信息 ----
    # init_guess='1e'：避免空基组 mol 下 'minao' 等 guess 报错
    # max_cycle=0  ：不做 SCF 迭代，只计算初始能量并填充必要属性
    mf.init_guess = "1e"
    mf.max_cycle = 0
    mf.kernel()

    # 部分 PySCF 版本在 max_cycle=0 时不会写入 mo_coeff，这里兜底。
    # FCIDUMP 的轨道本身已正交归一，单位矩阵是天然合理的初猜。
    if getattr(mf, "mo_coeff", None) is None:
        norb = mf.get_ovlp().shape[0]
        mf.mo_coeff = np.eye(norb)

    # ---- 3. 构建 MCSCF 对象 ----
    mc_type = mc_type.lower()
    if mc_type == "casscf":
        mc = mcscf.CASSCF(mf, ncas, nelecas)
    elif mc_type == "casci":
        mc = mcscf.CASCI(mf, ncas, nelecas)
    else:
        raise ValueError("mc_type must be 'casscf' or 'casci'")

    # ---- 4. 活性空间轨道选择 ----
    if active_orbs is not None:
        # sort_mo 将 mo_coeff 重排为 [core, active, virtual]
        # caslst 需包含 core + active 的索引
        mc.sort_mo(active_orbs)

    return mol, mf, mc
