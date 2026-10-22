# coding=UTF-8
"""用 F2M_iCIPT2 作为 CI solver 的 MCSCF driver。

F2M_iCIPT2 (``src/app/F2M_iCIPT2_Fe2S2.cpp``) 是 iCIPT2 家族中面向
[2Fe-2S] 模型 Hamiltonian 的 BUGA (Block Unitary Group Approach) 版本。
它读取

    <executable> <input> <FCIDUMP> [mx_spintwo] [mx_ext_level] [dm_flag]

其中 ``<input>`` 是与标准 iCIPT2 相同格式的输入文件（见
``iCIPT2_Loong/examples/06-F2M_iCIPT2/01-fe2s2/input``）。程序把收敛后
的 1-/2-RDM 写到当前目录下的

    fe2s2_dm1.csv , fe2s2_dm2.csv

（形如 ``state,i,j,rdm1`` / ``state,p,q,r,s,rdm2``，只存下三角）。

本模块提供两层接口：

1. :class:`F2MiCIPT2Solver`
   一个 PySCF FCI-solver 风格的 wrapper（``kernel`` / ``make_rdm1`` /
   ``make_rdm12`` / ``spin_square``），可以直接赋给
   ``pyscf.mcscf.CASSCF.fcisolver``。

2. :func:`kernel`
   MCSCF driver：从一个 FCIDUMP 出发构造 PySCF CASSCF 对象（借助
   :mod:`pyscf_util.misc.setup_mcscf_from_fcidump`），把 CI solver 换成
   F2M_iCIPT2，然后做宏迭代。

注意 / Note
------------
* ``F2M_iCIPT2_Fe2S2`` 内部把轨道顺序、``ref_occ``、``segment`` 全部
  写死为 20 轨道的 [2Fe-2S] 模型（轨道顺序见 ``01-fe2s2`` 例子）。
  因此本 driver 要求 **ncas == norb(FCIDUMP)**，并且保持 FCIDUMP 的
  轨道基（不做任何轨道旋转）。对 20 轨道 30 电子模型，这等价于
  full-valence CASSCF，轨道旋转是冗余的；driver 仍然走完整的 MCSCF
  流程（积分 dump -> CI -> RDM -> 梯度/Hessian），以保证接口正确。
* 由于轨道基不能旋转，这里关掉 PySCF 的 ``canonicalization``，否则
  宏迭代后 ``mo_coeff`` 会被重新酉旋转，F2M_iCIPT2 的局域轨道图像
  （Fe 3d / S 3p 分组）会被破坏。
* F2M_iCIPT2 的可执行文件通过环境变量 ``F2M_ICIPT2`` 指定；若未设置，
  会尝试在 ``PATH`` 中查找 ``F2M_iCIPT2_Fe2S2.exe``。
* F2M_iCIPT2 只在 ``argc >= 6`` 时计算 RDM，而此时 ``argv[4]`` 同时
  被当作 ``mx_ext_level`` 使用（程序里的小 bug）。所以 ``mx_ext_level
  > 0`` 才能拿到 RDM；``mx_ext_level == 0`` 时只能做能量计算。
"""

import os
import re
import shutil
import subprocess

import numpy

from pyscf import ao2mo, lib, tools
from pyscf.lib import logger

# F2M_iCIPT2_Fe2S2 写出的 RDM 文件名是写死的
F2M_DM1_NAME = "fe2s2_dm1.csv"
F2M_DM2_NAME = "fe2s2_dm2.csv"

# F2M_iCIPT2_Fe2S2 内部写死的模型：20 个轨道
F2M_NORB = 20

# 与 iCIPT2_Loong/examples/06-F2M_iCIPT2/01-fe2s2/input 保持一致
_DEFAULT_INPUT_TEMPLATE = """\
cmin={cmin}
etol={etol}
nsegment={nsegment}
nvalelec={nvalelec}
rotatemo=0
inputcfg=3
davidson=0
selection=0
dumprdm=2
relative=0
perturbation=0 0
task={task}
print=11
spin={spin}
nstates={nstates}
"""

# 解析 stdout 中的 CI 能量:  "Eig[   0]  :   -116.47094408"
_ENERGY_PATTERN = re.compile(
    r"Eig\[\s*(\d+)\s*\]\s*:\s*(-?\d+\.\d+(?:[eE][-+]?\d+)?)"
)


def _locate_executable():
    """找到 F2M_iCIPT2 可执行文件。"""
    exe = os.environ.get("F2M_ICIPT2")
    if exe:
        return exe
    exe = shutil.which("F2M_iCIPT2_Fe2S2.exe")
    if exe:
        return exe
    raise RuntimeError(
        "Cannot find the F2M_iCIPT2 executable. "
        "Please set the environment variable F2M_ICIPT2, e.g. "
        "export F2M_ICIPT2=/path/to/iCIPT2_Loong/bin/F2M_iCIPT2_Fe2S2.exe"
    )


###########################################################
# 积分 / RDM 文件的读写
###########################################################


def _write_fcidump(filename, h1e, eri, norb, nelec, ecore=0, ms=0, orbsym=None):
    """把活性空间积分写成 F2M_iCIPT2 能读的 FCIDUMP。

    ``h1e`` 是 PySCF 传进来的有效单电子积分（已经包含 core 贡献），
    ``ecore`` 是与之配套的常数项，二者一起写进 FCIDUMP。
    """
    eri8 = ao2mo.restore(8, eri, norb)
    tools.fcidump.from_integrals(
        filename,
        h1e,
        eri8,
        norb,
        nelec,
        ecore,
        ms=ms,
        orbsym=orbsym,
    )
    return filename


def _read_csv_blocks(filename):
    """读取 iCIPT2 风格的 RDM csv，返回 (has_state, rows)。

    rows 是 ``(state, i, j, ..., value)`` 形式的元组列表；表头有两种：

    * 不做态平均: ``state,i,j,rdm1`` / ``state,p,q,r,s,rdm2``
    * 做态平均  : ``i,j,rdm1``       / ``p,q,r,s,rdm2``
    """
    with open(filename) as f:
        header = f.readline().strip().replace(" ", "").split(",")
        has_state = header[0].lower() == "state"
        ncol = len(header)
        rows = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split(",")
            if len(cols) != ncol:
                raise ValueError(
                    "Bad line in %s: %s" % (filename, line)
                )
            rows.append(tuple(cols))
    return has_state, rows


def read_rdm1_csv(filename, norb, state=0):
    """读入 1-RDM，返回 (norb, norb) 的对称矩阵。"""
    has_state, rows = _read_csv_blocks(filename)
    rdm1 = numpy.zeros((norb, norb))
    for row in rows:
        if has_state:
            k, i, j, val = row
            if int(k) != state:
                continue
        else:
            i, j, val = row
        i, j, val = int(i), int(j), float(val)
        if i >= norb or j >= norb:
            raise ValueError("Orbital index out of range in %s" % filename)
        rdm1[i, j] = rdm1[j, i] = val
    return rdm1


def read_rdm2_csv(filename, norb, state=0):
    """读入 2-RDM 并转换成 PySCF 的存储顺序。

    iCIPT2 的 csv 只存 8 重对称意义下的非零元素，标签为 ``(p,q,r,s)``，
    数值对应 :math:`\\langle p^\\dagger q^\\dagger s r\\rangle`。先按
    :math:`\\Gamma_{pqrs}=\\Gamma_{qpsr}=\\Gamma_{rspq}=\\Gamma_{srqp}`
    补全，再做 ``transpose(0, 3, 2, 1)``，得到与
    :func:`pyscf.fci.direct_spin1.make_rdm12` 相同顺序的数组：

        dm2[p, q, r, s] -> 与 chemists' eri 收缩: 0.5 * einsum('pqrs,pqrs')

    该转换已用 "E = ecore + h1.RDM1 + 0.5 (pq|rs) RDM2" 以及
    "sum_r dm2[p,q,r,r] = (N-1) rdm1[p,q]" 两条恒等式验证。
    """
    has_state, rows = _read_csv_blocks(filename)

    p = numpy.empty(len(rows), dtype=int)
    q = numpy.empty(len(rows), dtype=int)
    r = numpy.empty(len(rows), dtype=int)
    s = numpy.empty(len(rows), dtype=int)
    val = numpy.empty(len(rows))

    nkeep = 0
    for row in rows:
        if has_state:
            k, pp, qq, rr, ss, vv = row
            if int(k) != state:
                continue
        else:
            pp, qq, rr, ss, vv = row
        p[nkeep] = int(pp)
        q[nkeep] = int(qq)
        r[nkeep] = int(rr)
        s[nkeep] = int(ss)
        val[nkeep] = float(vv)
        nkeep += 1

    p, q, r, s, val = p[:nkeep], q[:nkeep], r[:nkeep], s[:nkeep], val[:nkeep]
    if nkeep and (max(p.max(), q.max(), r.max(), s.max()) >= norb):
        raise ValueError("Orbital index out of range in %s" % filename)

    # 8 重对称补全
    rdm2 = numpy.zeros((norb, norb, norb, norb))
    rdm2[p, q, r, s] = val
    rdm2[q, p, s, r] = val
    rdm2[r, s, p, q] = val
    rdm2[s, r, q, p] = val

    # 转成 PySCF 的顺序
    return rdm2.transpose(0, 3, 2, 1)


def rdm_energy(h1e, eri, rdm1, rdm2, ecore=0.0):
    """用 RDM 计算能量: E = ecore + h1.RDM1 + 0.5 (pq|rs) RDM2。"""
    norb = h1e.shape[0]
    eri1 = ao2mo.restore(1, eri, norb)
    e1 = numpy.einsum("pq,pq->", h1e, rdm1)
    e2 = 0.5 * numpy.einsum("pqrs,pqrs->", eri1, rdm2)
    return ecore + e1 + e2


def parse_ci_energies(text, nstates):
    """从 F2M_iCIPT2 的 stdout 中取出最后一次对角化的各个态能量。"""
    matches = _ENERGY_PATTERN.findall(text)
    if len(matches) < nstates:
        raise RuntimeError(
            "Cannot find %d CI energies in the F2M_iCIPT2 output" % nstates
        )
    last = matches[-nstates:]
    energies = []
    for i, (istate, e) in enumerate(last):
        if int(istate) != i and nstates > 1:
            raise RuntimeError("Inconsistent state index in F2M_iCIPT2 output")
        energies.append(float(e))
    return energies


###########################################################
# FCI-solver 风格的 wrapper
###########################################################


class F2MiCIPT2Solver(lib.StreamObject):
    """把 F2M_iCIPT2 包装成 PySCF 的 FCI solver。

    Attributes
    ----------
    executable : str
        F2M_iCIPT2 可执行文件。
    workdir : str
        工作目录；FCIDUMP、input、stdout 和 RDM csv 都写在这里。
    spin : int
        2S（PySCF 约定），默认 0（singlet）。
    nroots : int
        态的个数。目前只支持 1，因为 F2M_iCIPT2 只 dump 一个态的 RDM。
    cmin : float or str
        iCIPT2 的 selection threshold。
    etol : float
        iCIPT2 的能量收敛阈值。
    mx_spintwo : int
        BUGA 局域态允许的最大自旋（程序 argv[3]）。
    mx_ext_level : int
        S 3p 外围配体局域激发等级（程序 argv[4]）。
    nsegment, nvalelec :
        与 iCIPT2 输入文件中的同名关键词一致（BUGA 版本里主要用于
        读取积分，默认值与官方 Fe2S2 例子相同）。
    energy_source : {'rdm', 'printed'}
        返回给 PySCF 的能量来源。'rdm'（默认）用 RDM 重算能量，保证
        能量和密度矩阵自洽；'printed' 用程序在 stdout 里最后打印的
        能量（final prune 之前）。正常计算两者一致。
    """

    def __init__(
        self,
        executable=None,
        workdir=None,
        spin=0,
        nroots=1,
        cmin=1e-4,
        etol=1e-10,
        mx_spintwo=10,
        mx_ext_level=1,
        nsegment="0 0 10 10 0 0",
        nvalelec=2,
        irrep=0,
        calc_rdm=True,
        energy_source="rdm",
        omp_num_threads=None,
        verbose=None,
    ):
        self.executable = executable if executable is not None else _locate_executable()
        self.workdir = os.path.abspath(workdir if workdir is not None else "f2m_icipt2_work")
        self.spin = spin
        self.nroots = nroots
        self.cmin = cmin
        self.etol = etol
        self.mx_spintwo = mx_spintwo
        self.mx_ext_level = mx_ext_level
        self.nsegment = nsegment
        self.nvalelec = nvalelec
        self.irrep = irrep
        self.calc_rdm = calc_rdm
        if energy_source not in ("rdm", "printed"):
            raise ValueError("energy_source must be 'rdm' or 'printed'")
        self.energy_source = energy_source
        self.omp_num_threads = omp_num_threads
        self.verbose = logger.NOTE if verbose is None else verbose

        self.converged = True
        self.e_ci = None
        self.e_printed = None
        self.mol = None

        self._rdm1 = None
        self._rdm2 = None
        self._cached_key = None
        self._cached_h1e = None
        self._cached_eri1 = None
        self._cached_result = None
        self._ncall = 0

        self._keys = set(self.__dict__.keys())

    #######################################################

    def dump_flags(self, verbose=None):
        log = logger.new_logger(self, verbose)
        log.info("")
        log.info("******** F2M_iCIPT2 solver ********")
        log.info("executable    = %s", self.executable)
        log.info("workdir       = %s", self.workdir)
        log.info("spin (2S)     = %d", self.spin)
        log.info("nroots        = %d", self.nroots)
        log.info("cmin / etol   = %s / %s", self.cmin, self.etol)
        log.info("mx_spintwo    = %d", self.mx_spintwo)
        log.info("mx_ext_level  = %d", self.mx_ext_level)
        log.info("nsegment      = %s", self.nsegment)
        log.info("energy source = %s", self.energy_source)
        log.info("OMP_NUM_THREADS = %s", self.omp_num_threads if self.omp_num_threads is not None else "(inherited)")
        log.info("")
        return self

    def _task_line(self, spintwo):
        return "%d %d %d %s" % (
            spintwo,
            self.irrep,
            self.nroots,
            " ".join(["1"] * self.nroots),
        )

    def _write_input(self, filename, spintwo):
        text = _DEFAULT_INPUT_TEMPLATE.format(
            cmin=self.cmin,
            etol=self.etol,
            nsegment=self.nsegment,
            nvalelec=self.nvalelec,
            task=self._task_line(spintwo),
            spin=spintwo / 2.0,
            nstates=self.nroots,
        )
        with open(filename, "w") as f:
            f.write(text)
        return filename

    def _run_program(self, inputfile, integralfile, spintwo):
        # argv: input, FCIDUMP, mx_spintwo, mx_ext_level[, dm_flag]
        # NOTE: F2M_iCIPT2 只在 argc >= 6 时算 RDM，且把 argv[4]
        #       (mx_ext_level) 当成 dm flag 读，见模块 docstring。
        cmd = [
            self.executable,
            inputfile,
            integralfile,
            str(self.mx_spintwo),
            str(self.mx_ext_level),
        ]
        if self.calc_rdm:
            if self.mx_ext_level <= 0:
                raise ValueError(
                    "F2M_iCIPT2 only dumps RDMs when mx_ext_level > 0 "
                    "(argv[4] is used as the RDM flag)"
                )
            cmd.append("1")

        env = os.environ.copy()
        if self.omp_num_threads is not None:
            env["OMP_NUM_THREADS"] = str(int(self.omp_num_threads))

        logfile = os.path.join(self.workdir, "f2m_icipt2.out")
        errfile = os.path.join(self.workdir, "f2m_icipt2.err")
        with open(logfile, "w") as out, open(errfile, "w") as err:
            proc = subprocess.run(
                cmd, cwd=self.workdir, env=env, stdout=out, stderr=err
            )

        with open(logfile) as f:
            output = f.read()

        if proc.returncode != 0:
            tail = "\n".join(output.splitlines()[-30:])
            raise RuntimeError(
                "F2M_iCIPT2 failed (return code %d).\n"
                "Last lines of %s:\n%s" % (proc.returncode, logfile, tail)
            )
        return output, logfile

    #######################################################

    def kernel(self, h1e, eri, norb, nelec, ci0=None, ecore=0, **kwargs):
        """运行 F2M_iCIPT2，返回 (e_ci + ecore, ci)。"""
        h1e = numpy.asarray(h1e)
        if isinstance(nelec, (int, numpy.integer)):
            nelecb = nelec // 2
            neleca = nelec - nelecb
        else:
            neleca, nelecb = nelec
        nelec_tot = neleca + nelecb
        spintwo = abs(neleca - nelecb)

        if self.nroots != 1:
            raise NotImplementedError(
                "F2MiCIPT2Solver only supports a single state (nroots=1)"
            )
        if norb != F2M_NORB:
            raise ValueError(
                "F2M_iCIPT2_Fe2S2 is compiled for the %d-orbital [2Fe-2S] "
                "model, but norb = %d was passed.  All %d fcidump orbitals "
                "must be active." % (F2M_NORB, norb, F2M_NORB)
            )
        if self.spin is not None and int(self.spin) != spintwo:
            raise ValueError(
                "Spin of the solver (%d) is inconsistent with nelec = %s"
                % (self.spin, str(nelec))
            )

        # 积分没变时直接复用上次结果（宏迭代里 PySCF 会用相同积分
        # 调用 solver 若干次）
        eri1 = numpy.asarray(ao2mo.restore(1, eri, norb))
        key = (h1e.shape, float(ecore), nelec_tot, spintwo)
        if (
            key == self._cached_key
            and self._cached_h1e is not None
            and numpy.allclose(h1e, self._cached_h1e, rtol=0, atol=1e-12)
            and numpy.allclose(eri1, self._cached_eri1, rtol=0, atol=1e-12)
        ):
            e_tot, ci = self._cached_result
            logger.debug(self, "F2M_iCIPT2: reuse the previous result")
            return e_tot, ci.copy()

        os.makedirs(self.workdir, exist_ok=True)
        for fname in (F2M_DM1_NAME, F2M_DM2_NAME):
            fpath = os.path.join(self.workdir, fname)
            if os.path.isfile(fpath):
                os.remove(fpath)

        integralfile = os.path.join(self.workdir, "FCIDUMP")
        _write_fcidump(
            integralfile,
            h1e,
            eri,
            norb,
            nelec_tot,
            ecore=ecore,
            ms=spintwo,
        )

        inputfile = os.path.join(self.workdir, "f2m_icipt2.inp")
        self._write_input(inputfile, spintwo)

        self._ncall += 1
        output, logfile = self._run_program(inputfile, integralfile, spintwo)

        energies = parse_ci_energies(output, self.nroots)
        # FCIDUMP 里已经带了 ecore，所以打印出来的能量就是总能量
        self.e_printed = energies[0]

        if not self.calc_rdm:
            dm1 = dm2 = None
            e_rdm = None
        else:
            dm1_path = os.path.join(self.workdir, F2M_DM1_NAME)
            dm2_path = os.path.join(self.workdir, F2M_DM2_NAME)
            if not (os.path.isfile(dm1_path) and os.path.isfile(dm2_path)):
                raise RuntimeError(
                    "F2M_iCIPT2 did not dump the RDMs in %s" % self.workdir
                )
            dm1 = read_rdm1_csv(dm1_path, norb)
            dm2 = read_rdm2_csv(dm2_path, norb)

            # sanity check 1: 电子数
            if abs(numpy.trace(dm1) - nelec_tot) > 1e-6:
                raise RuntimeError(
                    "Bad 1-RDM from F2M_iCIPT2: trace = %.6f, expected %d"
                    % (numpy.trace(dm1), nelec_tot)
                )
            # sanity check 2: 2-RDM 收缩回 1-RDM
            if nelec_tot > 1:
                dm1_from_2 = numpy.einsum("pqrr->pq", dm2) / (nelec_tot - 1)
                err = abs(dm1_from_2 - dm1).max()
                if err > 1e-6:
                    raise RuntimeError(
                        "Inconsistent 1-/2-RDM from F2M_iCIPT2 (err = %.3e)"
                        % err
                    )

            e_rdm = rdm_energy(h1e, eri, dm1, dm2, ecore)

            # F2M_iCIPT2 在 final prune 之后不再打印能量；当 cmin 很小、
            # 空间很大时，dump 出来的 RDM 所对应的态可能与最后打印的
            # 能量不一致（实测 cmin=1e-4 时约 0.3 eV）。为保证 MCSCF 的
            # 能量与 RDM 自洽，这里统一采用 RDM 能量，并给出 warning。
            if abs(e_rdm - self.e_printed) > 1e-6:
                logger.warn(
                    self,
                    "F2M_iCIPT2: energy printed by the program (%.12f) "
                    "differs from the RDM energy (%.12f) by %.3e. "
                    "The RDM energy will be used.",
                    self.e_printed,
                    e_rdm,
                    e_rdm - self.e_printed,
                )

        if self.calc_rdm and self.energy_source == "rdm":
            e_ci = e_rdm
        else:
            e_ci = self.e_printed
        self.e_ci = e_ci
        self._rdm1, self._rdm2 = dm1, dm2

        # PySCF 用 fcivec 作为 solvers 之间的句柄；这里没有 CI 向量，
        # 返回一个长度 1 的数组作为占位符。
        ci = numpy.array([e_ci])

        self._cached_key = key
        self._cached_h1e = h1e.copy()
        self._cached_eri1 = eri1.copy()
        self._cached_result = (e_ci, ci.copy())
        self.converged = True

        logger.info(
            self,
            "F2M_iCIPT2 call %d: E(printed) = %.12f  E(rdm) = %s  log = %s",
            self._ncall,
            self.e_printed,
            "%.12f" % e_rdm if e_rdm is not None else "N/A",
            logfile,
        )
        return e_ci, ci

    def make_rdm1(self, fcivec=None, norb=None, nelec=None, **kwargs):
        if self._rdm1 is None:
            raise RuntimeError(
                "No RDM available.  Set calc_rdm=True and call kernel() first."
            )
        return self._rdm1

    def make_rdm2(self, fcivec=None, norb=None, nelec=None, **kwargs):
        if self._rdm2 is None:
            raise RuntimeError(
                "No RDM available.  Set calc_rdm=True and call kernel() first."
            )
        return self._rdm2

    def make_rdm12(self, fcivec=None, norb=None, nelec=None, **kwargs):
        return self.make_rdm1(), self.make_rdm2()

    def spin_square(self, fcivec=None, norb=None, nelec=None):
        s = float(self.spin) / 2.0
        return s * (s + 1.0), self.spin + 1


###########################################################
# MCSCF driver
###########################################################


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
    mx_spintwo=10,
    mx_ext_level=1,
    nsegment="0 0 10 10 0 0",
    nvalelec=2,
    calc_rdm=True,
    energy_source="rdm",
    omp_num_threads=None,
    mc_conv_tol=1e-7,
    mc_max_macro=32,
    run_mcscf=True,
    verbose=None,
):
    """MCSCF driver: FCIDUMP -> PySCF CASSCF，CI solver 用 F2M_iCIPT2。

    Parameters
    ----------
    fcidump : str
        模型 Hamiltonian 的 FCIDUMP 文件（例如
        ``Active-space-model-for-Iron-Sulfur-Clusters/Fe2S2_and_Fe4S4/Fe2S2/fe2s2``）。
    ncas, nelecas : int, optional
        活性轨道数和活性电子数。默认取 FCIDUMP 里的全部轨道 / 电子。
        F2M_iCIPT2_Fe2S2 只能处理 20 轨道模型，所以 ``ncas`` 必须等于
        FCIDUMP 的 NORB。
    active_orbs : None
        为了保持 F2M_iCIPT2 写死的轨道顺序，这里不允许重排轨道。
    spin : int
        2S，默认 0。
    cmin, etol, mx_spintwo, mx_ext_level, nsegment, nvalelec :
        透传给 F2M_iCIPT2（见 :class:`F2MiCIPT2Solver`）。
    energy_source : {'rdm', 'printed'}
        透传给 F2M_iCIPT2（见 :class:`F2MiCIPT2Solver`）。
    omp_num_threads : int, optional
        给 F2M_iCIPT2 子进程设置的 OMP_NUM_THREADS；None 表示继承环境。
    run_mcscf : bool
        True  : 走完整的 CASSCF 宏迭代；
        False : 只做一次 CASCI（一次 F2M_iCIPT2 调用）。

    Returns
    -------
    mc : pyscf.mcscf.CASSCF
    """
    from pyscf_util.misc.setup_mcscf_from_fcidump import setup_mcscf_from_fcidump

    data = tools.fcidump.read(fcidump)
    norb = data["NORB"]
    nelec_fcidump = data["NELEC"]

    if ncas is None:
        ncas = norb
    if nelecas is None:
        nelecas = nelec_fcidump

    if ncas != norb:
        raise ValueError(
            "F2M_iCIPT2_Fe2S2 requires all %d fcidump orbitals to be active "
            "(ncas = %d was requested).  Its local-orbital structure and "
            "ref_occ are hard-coded for the full model."
            % (norb, ncas)
        )
    if active_orbs is not None:
        raise ValueError(
            "active_orbs is not supported: F2M_iCIPT2 requires the original "
            "fcidump orbital order."
        )
    if norb != F2M_NORB:
        raise ValueError(
            "F2M_iCIPT2_Fe2S2 is compiled for %d orbitals, this fcidump has %d"
            % (F2M_NORB, norb)
        )
    if run_mcscf and not calc_rdm:
        raise ValueError(
            "MCSCF needs the RDMs of the CI solver; set calc_rdm=True."
        )

    log = logger.new_logger(None, verbose)
    log.info("")
    log.info("******** MCSCF with F2M_iCIPT2 solver ********")
    log.info("fcidump = %s", fcidump)
    log.info("norb    = %d, nel = %d, ncas = %d, nelecas = %s",
             norb, nelec_fcidump, ncas, str(nelecas))

    mol, mf, mc = setup_mcscf_from_fcidump(
        fcidump,
        ncas=ncas,
        nelecas=nelecas,
        active_orbs=None,
        mc_type="casscf",
        verbose=0,
    )

    # 关键: 保持 FCIDUMP 的轨道基。FCIDUMP 的轨道本身正交，
    # 单位矩阵就是最自然的 mo_coeff。
    mo_init = numpy.eye(norb)
    mc.mo_coeff = mo_init

    solver = F2MiCIPT2Solver(
        executable=executable,
        workdir=workdir,
        spin=spin,
        nroots=1,
        cmin=cmin,
        etol=etol,
        mx_spintwo=mx_spintwo,
        mx_ext_level=mx_ext_level,
        nsegment=nsegment,
        nvalelec=nvalelec,
        calc_rdm=calc_rdm,
        energy_source=energy_source,
        omp_num_threads=omp_num_threads,
        verbose=verbose,
    )
    mc.fcisolver = solver

    mc.conv_tol = mc_conv_tol
    mc.max_cycle_macro = mc_max_macro
    # 不做轨道酉旋转，否则 F2M_iCIPT2 的局域轨道图像会被破坏
    mc.canonicalization = False

    if run_mcscf:
        mc.kernel(mo_init)
        mo_final = mc.mo_coeff
        # mo_coeff 必须仍然是 FCIDUMP 轨道（允许行置换/符号，不允许混合）
        offdiag = numpy.abs(mo_final) - numpy.abs(numpy.eye(norb))
        if offdiag.max() > 1e-6:
            logger.warn(
                mc,
                "mo_coeff is no longer the fcidump orbital basis "
                "(max |offdiag| = %.3e).  The F2M_iCIPT2 result is only "
                "meaningful in the original local orbital basis.",
                offdiag.max(),
            )
    else:
        mc.e_tot, mc.e_cas, mc.ci = mc.casci(mo_init)

    log.info("")
    log.info("MCSCF/F2M_iCIPT2 total energy = %.12f", mc.e_tot)
    log.info("MCSCF/F2M_iCIPT2 CI    energy = %.12f", mc.e_cas)
    return mc


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        fcidump = sys.argv[1]
    else:
        # pyscf_util/MeanField/ -> pyscf_util2/ -> GithubRepo/
        _repo = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "..", "..")
        )
        fcidump = os.path.join(
            _repo,
            "Active-space-model-for-Iron-Sulfur-Clusters",
            "Fe2S2_and_Fe4S4",
            "Fe2S2",
            "fe2s2",
        )
        if not os.path.isfile(fcidump):
            raise SystemExit(
                "FCIDUMP not found: %s\n\nUsage: python f2micipt2.py <FCIDUMP>"
                % fcidump
            )

    mc = kernel(fcidump, cmin=1e-2, etol=1e-6, verbose=4)
    print("E = ", mc.e_tot)
