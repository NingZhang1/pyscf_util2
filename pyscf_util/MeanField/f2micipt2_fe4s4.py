# coding=UTF-8
"""用 ``F2M_iCIPT2_Fe4S4`` 作为 CI solver 的 MCSCF driver。

``F2M_iCIPT2_Fe4S4`` 是 ``F2M_iCIPT2_Fe2S2`` 的 [4Fe-4S] 版本
(``src/app/F2M_iCIPT2_Fe4S4.cpp``)：36 个轨道、54 个电子，轨道顺序为

    0, 1                     S  (peripheral)
    2 .. 6                   Fe 1
    7 .. 11                  Fe 2
    12 .. 23                 S  (4 x 3p, bridging)
    24 .. 28                 Fe 3
    29 .. 33                 Fe 4
    34, 35                   S  (peripheral)

局域态模型支持三种累加模式：

============  ==========================================================
mode 1        spin excitation only: Fe(3+) (d^5) / Fe(2+) (d^6)
mode 2        + MMCT: 额外 Fe(+) (d^7) / Fe(4+) (d^4, high spin)
mode 3        + LMCT: 配体 S 3p 的局域激发，激发等级和自旋上限由
              ``mx_ligand_ext_level`` / ``mx_ligand_spintwo`` 控制
============  ==========================================================

程序接口::

    <executable> <input> <FCIDUMP> <mx_spintwo> <mode> \\
                 <mx_ligand_ext_level> <mx_ligand_spintwo> <calculate_dm>

其中 ``<input>`` 是与 iCIPT2 相同格式的输入文件。收敛后的 1-/2-RDM
写到当前目录的 ``fe4s4_dm1.csv`` / ``fe4s4_dm2.csv``。

本模块只是 :mod:`pyscf_util.MeanField.f2micipt2` 的一个子类：模型相关的
常量（轨道数、RDM 文件名、默认 ``nsegment`` / ``mx_spintwo``、可执行文件
环境变量）通过类属性覆盖。用法见
``examples/15-F2M-iCIPT2-Fe4S4-MCSCF.py``。
"""

from pyscf_util.MeanField import f2micipt2 as _fe2s2
from pyscf_util.MeanField.f2micipt2 import (
    _write_fcidump,  # noqa: F401  (re-export for convenience)
    extract_important_csf,  # noqa: F401
    parse_ci_energies,  # noqa: F401
    parse_local_states,  # noqa: F401
    rdm_energy,  # noqa: F401
    read_rdm1_csv,  # noqa: F401
    read_rdm2_csv,  # noqa: F401
)

# F2M_iCIPT2_Fe4S4 写出的 RDM 文件名是写死的
F2M_DM1_NAME = "fe4s4_dm1.csv"
F2M_DM2_NAME = "fe4s4_dm2.csv"

# F2M_iCIPT2_Fe4S4 内部写死的模型：36 个轨道
F2M_NORB = 36

# 三种局域态模式
MODES = {
    1: "spin excitation only (Fe2+/Fe3+)",
    2: "MMCT (Fe+/Fe4+ included)",
    3: "LMCT (ligand excitations included)",
}


class F2MiCIPT2Fe4S4Solver(_fe2s2.F2MiCIPT2Solver):
    """把 ``F2M_iCIPT2_Fe4S4`` 包装成 PySCF 的 FCI solver。

    除了基类的参数外，``mode`` / ``mx_ligand_ext_level`` /
    ``mx_ligand_spintwo`` 的含义见模块 docstring。
    """

    model_name = "Fe4S4"
    model_tag = "Fe4S4"
    norb_model = F2M_NORB
    dm1_name = F2M_DM1_NAME
    dm2_name = F2M_DM2_NAME
    default_nsegment = "0 0 18 18 0 0"
    default_mx_spintwo = 20
    exe_env = "F2M_ICIPT2_FE4S4"
    exe_name = "F2M_iCIPT2_Fe4S4.exe"


def kernel(
    fcidump,
    ncas=None,
    nelecas=None,
    active_orbs=None,
    executable=None,
    workdir=None,
    spin=0,
    cmin=1e-4,
    etol=1e-10,
    mx_spintwo=None,
    mode=1,
    mx_ligand_ext_level=1,
    mx_ligand_spintwo=1,
    nsegment=None,
    nvalelec=2,
    calc_rdm=True,
    energy_source="rdm",
    omp_num_threads=None,
    print_ci_each_call=True,
    n_top_csf=3,
    mc_conv_tol=1e-7,
    mc_max_macro=32,
    run_mcscf=True,
    internal_rotation=True,
    print_summary=False,
    dump_mo=True,
    solver=None,
    verbose=None,
):
    """MCSCF driver: Fe4S4 FCIDUMP -> PySCF CASSCF，CI solver 用 F2M_iCIPT2_Fe4S4。

    参数含义与 :func:`pyscf_util.MeanField.f2micipt2.kernel` 相同，另外：

    mode : int
        1: 只有自旋激发 (Fe2+/Fe3+)；
        2: 允许 MMCT (Fe+/Fe4+)；
        3: 允许 LMCT (配体激发，自旋上限由 ``mx_ligand_spintwo`` 控制)。
    mx_ligand_ext_level, mx_ligand_spintwo :
        mode 3 下配体局域激发的等级和最大局域自旋。

    Returns
    -------
    mc : pyscf.mcscf.CASSCF
    """
    if solver is None:
        solver = F2MiCIPT2Fe4S4Solver(
            executable=executable,
            workdir=workdir,
            spin=spin,
            nroots=1,
            cmin=cmin,
            etol=etol,
            mx_spintwo=mx_spintwo,
            mode=mode,
            mx_ligand_ext_level=mx_ligand_ext_level,
            mx_ligand_spintwo=mx_ligand_spintwo,
            nsegment=nsegment,
            nvalelec=nvalelec,
            calc_rdm=calc_rdm,
            energy_source=energy_source,
            omp_num_threads=omp_num_threads,
            print_ci_each_call=print_ci_each_call,
            n_top_csf=n_top_csf,
            verbose=verbose,
        )

    return _fe2s2.kernel(
        fcidump,
        ncas=ncas,
        nelecas=nelecas,
        active_orbs=active_orbs,
        executable=executable,
        workdir=workdir,
        spin=spin,
        cmin=cmin,
        etol=etol,
        mx_spintwo=mx_spintwo,
        mode=mode,
        mx_ligand_ext_level=mx_ligand_ext_level,
        mx_ligand_spintwo=mx_ligand_spintwo,
        nsegment=nsegment,
        nvalelec=nvalelec,
        calc_rdm=calc_rdm,
        energy_source=energy_source,
        omp_num_threads=omp_num_threads,
        print_ci_each_call=print_ci_each_call,
        n_top_csf=n_top_csf,
        mc_conv_tol=mc_conv_tol,
        mc_max_macro=mc_max_macro,
        run_mcscf=run_mcscf,
        internal_rotation=internal_rotation,
        print_summary=print_summary,
        dump_mo=dump_mo,
        solver=solver,
        verbose=verbose,
    )


if __name__ == "__main__":
    import os
    import sys

    if len(sys.argv) > 1:
        fcidump = sys.argv[1]
    else:
        _repo = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "..", "..")
        )
        fcidump = os.path.join(
            _repo,
            "Active-space-model-for-Iron-Sulfur-Clusters",
            "Fe2S2_and_Fe4S4",
            "Fe4S4",
            "fe4s4",
        )
        if not os.path.isfile(fcidump):
            raise SystemExit(
                "FCIDUMP not found: %s\n\nUsage: python f2micipt2_fe4s4.py <FCIDUMP>"
                % fcidump
            )

    mc = kernel(fcidump, cmin=1e-1, etol=1e-4, mode=1, verbose=4,
                print_summary=True)
    print("E = ", mc.e_tot)
