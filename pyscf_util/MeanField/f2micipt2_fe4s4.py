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
mode 3        + LMCT: 配体 S 3p 的局域激发，激发态自旋上限由
              ``mx_ligand_spintwo`` 控制
============  ==========================================================

程序接口::

    <executable> <input> <FCIDUMP> <mx_spintwo> <mode> \\
                 <mx_ligand_ext_level> <mx_ligand_spintwo> <calculate_dm>

其中 ``<input>`` 是与 iCIPT2 相同格式的输入文件。收敛后的 1-/2-RDM
写到当前目录的 ``fe4s4_dm1.csv`` / ``fe4s4_dm2.csv``。

与 Fe2S2 版本一样，``F2M_iCIPT2_Fe4S4`` 的局域轨道顺序、``ref_occ``
都是写死的，所以这里要求 ``ncas == norb(FCIDUMP)`` 且保持 FCIDUMP 的
轨道基（不做轨道旋转，关闭 ``canonicalization``）。

用法见 ``examples/15-F2M-iCIPT2-Fe4S4-MCSCF.py``。
"""

import os
import shutil
import subprocess

from pyscf.lib import logger

from pyscf_util.MeanField import f2micipt2 as _fe2s2
from pyscf_util.MeanField.f2micipt2 import (
    _write_fcidump,
    parse_ci_energies,
    rdm_energy,
    read_rdm1_csv,
    read_rdm2_csv,
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


def _locate_executable():
    """找到 F2M_iCIPT2_Fe4S4 可执行文件。"""
    exe = os.environ.get("F2M_ICIPT2_FE4S4")
    if exe:
        return exe
    exe = os.environ.get("F2M_ICIPT2")
    if exe and os.path.basename(exe).startswith("F2M_iCIPT2_Fe4S4"):
        return exe
    exe = shutil.which("F2M_iCIPT2_Fe4S4.exe")
    if exe:
        return exe
    raise RuntimeError(
        "Cannot find the F2M_iCIPT2_Fe4S4 executable. "
        "Please set the environment variable F2M_ICIPT2_FE4S4, e.g. "
        "export F2M_ICIPT2_FE4S4=/path/to/iCIPT2_Loong/bin/F2M_iCIPT2_Fe4S4.exe"
    )


class F2MiCIPT2Fe4S4Solver(_fe2s2.F2MiCIPT2Solver):
    """把 ``F2M_iCIPT2_Fe4S4`` 包装成 PySCF 的 FCI solver。

    在 Fe2S2 版本的基础上增加了

    * ``mode``                : 1 / 2 / 3（见模块 docstring）
    * ``mx_ligand_ext_level`` : 配体局域激发等级 (mode 3 使用)
    * ``mx_ligand_spintwo``   : 配体激发态允许的最大局域 spintwo
    """

    def __init__(
        self,
        executable=None,
        workdir=None,
        spin=0,
        nroots=1,
        cmin=1e-4,
        etol=1e-10,
        mx_spintwo=20,
        mode=1,
        mx_ligand_ext_level=1,
        mx_ligand_spintwo=1,
        nsegment="0 0 18 18 0 0",
        nvalelec=2,
        irrep=0,
        calc_rdm=True,
        energy_source="rdm",
        omp_num_threads=None,
        verbose=None,
    ):
        if mode not in MODES:
            raise ValueError("mode must be 1, 2 or 3 (got %r)" % (mode,))
        if executable is None:
            executable = _locate_executable()
        super().__init__(
            executable=executable,
            workdir=workdir,
            spin=spin,
            nroots=nroots,
            cmin=cmin,
            etol=etol,
            mx_spintwo=mx_spintwo,
            mx_ext_level=mx_ligand_ext_level,
            nsegment=nsegment,
            nvalelec=nvalelec,
            irrep=irrep,
            calc_rdm=calc_rdm,
            energy_source=energy_source,
            omp_num_threads=omp_num_threads,
            verbose=verbose,
        )
        self.norb_model = F2M_NORB
        self.dm1_name = F2M_DM1_NAME
        self.dm2_name = F2M_DM2_NAME

        self.mode = mode
        self.mx_ligand_ext_level = mx_ligand_ext_level
        self.mx_ligand_spintwo = mx_ligand_spintwo

    def dump_flags(self, verbose=None):
        log = logger.new_logger(self, verbose)
        log.info("")
        log.info("******** F2M_iCIPT2_Fe4S4 solver ********")
        log.info("executable         = %s", self.executable)
        log.info("workdir            = %s", self.workdir)
        log.info("spin (2S)          = %d", self.spin)
        log.info("nroots             = %d", self.nroots)
        log.info("cmin / etol        = %s / %s", self.cmin, self.etol)
        log.info("mx_spintwo         = %d", self.mx_spintwo)
        log.info("mode               = %d (%s)", self.mode, MODES[self.mode])
        log.info("mx_ligand_ext      = %d", self.mx_ligand_ext_level)
        log.info("mx_ligand_spintwo  = %d", self.mx_ligand_spintwo)
        log.info("nsegment           = %s", self.nsegment)
        log.info("energy source      = %s", self.energy_source)
        log.info("")
        return self

    def _build_cmd(self, inputfile, integralfile):
        """F2M_iCIPT2_Fe4S4 的命令行参数。"""
        return [
            self.executable,
            inputfile,
            integralfile,
            str(self.mx_spintwo),
            str(self.mode),
            str(self.mx_ligand_ext_level),
            str(self.mx_ligand_spintwo),
            "1" if self.calc_rdm else "0",
        ]

    def _run_program(self, inputfile, integralfile, spintwo):
        cmd = self._build_cmd(inputfile, integralfile)

        env = os.environ.copy()
        if self.omp_num_threads is not None:
            env["OMP_NUM_THREADS"] = str(int(self.omp_num_threads))

        logfile = os.path.join(self.workdir, "f2m_icipt2_fe4s4.out")
        errfile = os.path.join(self.workdir, "f2m_icipt2_fe4s4.err")
        with open(logfile, "w") as out, open(errfile, "w") as err:
            proc = subprocess.run(
                cmd, cwd=self.workdir, env=env, stdout=out, stderr=err
            )

        with open(logfile) as f:
            output = f.read()

        if proc.returncode != 0:
            tail = "\n".join(output.splitlines()[-30:])
            raise RuntimeError(
                "F2M_iCIPT2_Fe4S4 failed (return code %d).\n"
                "Last lines of %s:\n%s" % (proc.returncode, logfile, tail)
            )
        return output, logfile


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
    mx_spintwo=20,
    mode=1,
    mx_ligand_ext_level=1,
    mx_ligand_spintwo=1,
    nsegment="0 0 18 18 0 0",
    nvalelec=2,
    calc_rdm=True,
    energy_source="rdm",
    omp_num_threads=None,
    mc_conv_tol=1e-7,
    mc_max_macro=32,
    run_mcscf=True,
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
        mx_ext_level=mx_ligand_ext_level,
        nsegment=nsegment,
        nvalelec=nvalelec,
        calc_rdm=calc_rdm,
        energy_source=energy_source,
        omp_num_threads=omp_num_threads,
        mc_conv_tol=mc_conv_tol,
        mc_max_macro=mc_max_macro,
        run_mcscf=run_mcscf,
        solver=solver,
        verbose=verbose,
    )


if __name__ == "__main__":
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

    mc = kernel(fcidump, cmin=1e-1, etol=1e-4, mode=1, verbose=4)
    print("E = ", mc.e_tot)
