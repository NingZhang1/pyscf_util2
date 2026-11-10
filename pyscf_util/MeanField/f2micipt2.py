# coding=UTF-8
"""用 F2M_iCIPT2 作为 CI solver 的 MCSCF driver。

``F2M_iCIPT2`` 是 iCIPT2 家族中面向 [2Fe-2S] / [4Fe-4S] 模型 Hamiltonian 的
BUGA (Block Unitary Group Approach) 版本 (``src/app/F2M_iCIPT2_Fe2S2.cpp`` /
``F2M_iCIPT2_Fe4S4.cpp``)。当前版本的程序接口是::

    <executable> <input> <FCIDUMP> <mx_spintwo> <mode> \\
                 <mx_ligand_ext_level> <mx_ligand_spintwo> <calculate_dm> \\
                 <only_ferro_antiferro> <hmat_brute_force_threshold> \\
                 <hmat_update> <prune_init_guess>

其中

* ``mx_spintwo``                : 允许的最大中间自旋 (Fe2S2 默认 10, Fe4S4 默认 20);
* ``mode``                      : 1 = 只有自旋激发 (Fe2+/Fe3+)；
                                  2 = + MMCT (Fe+/Fe4+ d^4 high spin)；
                                  3 = + LMCT (配体 S 3p 局域激发)；
* ``mx_ligand_ext_level``       : 配体局域激发等级 (mode 3)，S3p 支持 0..6；
* ``mx_ligand_spintwo``         : 配体激发态允许的最大 2S；
* ``calculate_dm``              : 已废弃（程序忽略），RDM 输出改由输入文件的
                                  ``dumprdm`` 控制；
* ``only_ferro_antiferro``      : 1 = 只保留铁磁(最大自旋)/反铁磁(最小自旋)耦合；
* ``hmat_brute_force_threshold``: > 0 时，维度小于该值时用 brute-force H；
* ``hmat_update``               : 1 = 增量式建 H (append-only)，0 = 每次重建；
* ``prune_init_guess``          : 1 = selection 时允许剪掉初猜 CSF。

输入文件与标准 iCIPT2 相同（见 ``iCIPT2_Loong/examples/06-F2M_iCIPT2/``），
其中两个关键字现在很重要：

* ``selection = 0`` : 不做 selection，初猜就是最终结果；
  ``selection = 1`` : 做 selection（由 ``cmin`` 控制）；
* ``dumprdm = 0/1/2/7`` : 0 = 不 dump；1 = 态平均 1RDM；
  2 = 态平均 1+2RDM；7 = 每个态分别 dump 1+2RDM。

收敛后的 1-/2-RDM 写到工作目录（Fe2S2: ``fe2s2_dm*.csv``,
Fe4S4: ``fe4s4_dm*.csv``）。

本模块提供两层接口：

1. :class:`F2MiCIPT2Solver`
   一个 PySCF FCI-solver 风格的 wrapper（``kernel`` / ``make_rdm1`` /
   ``make_rdm12`` / ``spin_square`` / ``print_summary`` / ``dump_mo``）。
   Fe4S4 版本见 :mod:`pyscf_util.MeanField.f2micipt2_fe4s4`。

2. :func:`kernel`
   MCSCF driver：从一个 FCIDUMP 出发构造 PySCF CASSCF 对象，把 CI solver
   换成 F2M_iCIPT2，然后做宏迭代。

注意 / Note
-----------
* 程序内部把轨道顺序、``ref_occ`` 写死，因此要求
  **ncas == norb(FCIDUMP)**，初猜使用 FCIDUMP 的轨道基。
* 默认打开 ``internal_rotation``，允许活性空间内部的轨道旋转，
  这样宏迭代会真正进行。F2M_iCIPT2 的局域态模型是在给定轨道基下定义的，
  旋转后的能量会随基变化；若只想在原始基下做 CASCI，设
  ``run_mcscf=False`` 或 ``internal_rotation=False``。
* 每次真正调用 F2M_iCIPT2 都会打印最重要的几个组态（``n_top_csf``），
  ``print_summary=True`` 时结束打印宏迭代历史、最后一次调用的局域态/
  重要组态、能量和 ``mo_coeff``；``dump_mo=True``（默认）把 ``mo_coeff``
  写成 ``mo_coeff.csv`` 并画出热图 ``mo_coeff_heatmap.png``。
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
selection={selection}
dumprdm={dumprdm}
relative=0
perturbation=0 0
task={task}
print=11
spin={spin}
nstates={nstates}
"""

# 解析 stdout 中的 CI 能量:  "Eig[   0]  :   -116.47094408"
_ENERGY_PATTERN = re.compile(r"Eig\[\s*(\d+)\s*\]\s*:\s*(-?\d+\.\d+(?:[eE][-+]?\d+)?)")

# 新版本程序的最终能量:  "task 0 final energies:" 后面跟
#   "  state  0 : -115.953720058215"
_TASK_HEADER_RE = re.compile(r"^task\s+(\d+)\s+final energies:\s*$", re.M)
_STATE_ENERGY_RE = re.compile(
    r"^\s*state\s+(\d+)\s*:\s*(-?\d+\.\d+(?:[eE][-+]?\d+)?)\s*$", re.M
)

_NCSF_RE = re.compile(r"ncsf\s*=\s*(\d+)")
_BLOCK_STATE_RE = re.compile(
    r"nelec : (\d+), spintwo : (\d+), nopen : (\d+)\s*\n\s*n_local_state\s*:\s*(\d+)"
)


def _locate_executable(
    exe_env="F2M_ICIPT2_FE2S2",
    exe_name="F2M_iCIPT2_Fe2S2.exe",
    model_tag="Fe2S2",
):
    """找到 F2M_iCIPT2 可执行文件。"""
    for var in (exe_env, "F2M_ICIPT2"):
        exe = os.environ.get(var)
        if exe and model_tag in os.path.basename(exe):
            return exe
    exe = os.environ.get(exe_env)
    if exe:
        return exe
    exe = shutil.which(exe_name)
    if exe:
        return exe
    raise RuntimeError(
        "Cannot find the F2M_iCIPT2 executable for %s. "
        "Please set the environment variable %s (or F2M_ICIPT2), e.g. "
        "export %s=/path/to/iCIPT2_Loong/bin/%s"
        % (model_tag, exe_env, exe_env, exe_name)
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
                raise ValueError("Bad line in %s: %s" % (filename, line))
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


###########################################################
# stdout 解析
###########################################################


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


def parse_task_energies(text, nstates=None):
    """解析 ``task N final energies:`` 段（新版本程序的最终能量）。

    返回最后一个 task 的各态能量；找不到时返回 ``None``。
    """
    headers = list(_TASK_HEADER_RE.finditer(text))
    if not headers:
        return None
    body = text[headers[-1].end() :]
    energies = [float(e) for _, e in _STATE_ENERGY_RE.findall(body)]
    if not energies:
        return None
    if nstates is not None:
        energies = energies[:nstates]
    return energies


def parse_local_states(text):
    """从 F2M_iCIPT2 的 stdout 中解析每个 block 的局域态列表。

    返回 ``[{'block': int, 'norb': int, 'states': [(nelec, 2S, nopen, nstate), ...]}, ...]``
    """
    blocks = []
    for chunk in text.split("PrimitiveLocalStateLoopUpTable")[1:]:
        m_block = re.search(r"block_id\s*:\s*(\d+)", chunk)
        m_norb = re.search(r"^\s*norb\s*:\s*(\d+)\s*$", chunk, re.M)
        if m_block is None:
            continue
        states = [
            (int(a), int(b), int(c), int(d))
            for a, b, c, d in _BLOCK_STATE_RE.findall(chunk)
        ]
        blocks.append(
            dict(
                block=int(m_block.group(1)),
                norb=int(m_norb.group(1)) if m_norb else None,
                states=states,
            )
        )
    return blocks


def parse_important_csfs(text):
    """解析最后一次 ``Important block csf`` 列表。

    返回 ``[{'label': str, 'coef': float, 'weight': float}, ...]``，
    weight 是日志里给出的百分比。
    """
    key = "Important block csf:"
    start = text.rfind(key)
    if start < 0:
        return []
    csfs = []
    for line in text[start:].splitlines()[1:]:
        line = line.rstrip()
        if not line or line.startswith("check importance") or line.startswith("*"):
            break
        m = re.search(r"(-?\d+\.\d+)\s*\(\s*([\d.]+)%\)\s*$", line)
        if m is None:
            continue
        label = line[: m.start()].strip()
        csfs.append(dict(label=label, coef=float(m.group(1)), weight=float(m.group(2))))
    return csfs


def extract_important_csf(text):
    """取出最后一次 ``PRINT BLOCK CSF SPACE`` 段（重要组态）。"""
    key = "Begin PRINT BLOCK CSF SPACE"
    start = text.rfind(key)
    if start < 0:
        return None
    end = text.find("End PRINT BLOCK CSF SPACE", start)
    if end < 0:
        return None
    # 去掉 "*************** End ..." 这一行
    end = text.rfind("\n", start, end) + 1
    return text[start:end].strip("\n")


def plot_matrix_heatmap(
    matrix,
    filename,
    title="",
    xlabel="MO index",
    ylabel="basis index",
    cmap="RdBu_r",
):
    """把矩阵画成热图并保存为图片 (需要 matplotlib)。"""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    matrix = numpy.asarray(matrix, dtype=float)
    nrow, ncol = matrix.shape
    vmax = float(numpy.abs(matrix).max())
    if vmax <= 0:
        vmax = 1.0

    fig, ax = plt.subplots(
        figsize=(max(6.0, 0.30 * ncol + 2.0), max(4.5, 0.30 * nrow + 1.5))
    )
    im = ax.imshow(
        matrix, cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto", origin="upper"
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xticks(range(0, ncol, max(1, ncol // 12)))
    ax.set_yticks(range(0, nrow, max(1, nrow // 12)))
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("coefficient")
    fig.tight_layout()
    fig.savefig(filename, dpi=150)
    plt.close(fig)
    return filename


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
        态的个数。目前只支持 1。
    cmin : float or str
        selection threshold。
    etol : float
        对角化/能量收敛阈值。
    mx_spintwo : int
        BUGA 局域态允许的最大中间自旋（程序 argv[3]）。
    mode : int
        1: 只有自旋激发 (Fe2+/Fe3+)；2: + MMCT (Fe+/Fe4+)；
        3: + LMCT (配体激发)。
    mx_ligand_ext_level, mx_ligand_spintwo :
        配体局域激发等级和最大局域自旋（mode 3；程序 argv[5], argv[6]）。
    selection : int
        输入文件关键字：0 = 不做 selection；1 = 做 selection。
    calc_rdm : bool
        是否请求 RDM。实际写什么由 ``dumprdm`` 决定：
        0 = 不写；1 = 态平均 1RDM；2 = 态平均 1+2RDM；7 = 每个态 1+2RDM。
    only_ferro_antiferro : bool
        只保留铁磁(最大自旋)/反铁磁(最小自旋)耦合。
    hmat_brute_force_threshold : int
        > 0 时，维度小于该值时用 brute-force 建 H，否则用 tree。
    hmat_update : bool
        True = 增量式建 H；False = 每次重建。
    prune_init_guess : bool
        selection 时是否允许剪掉初猜 CSF。
    energy_source : {'rdm', 'printed'}
        返回给 PySCF 的能量来源。'rdm'（默认）用 RDM 重算能量，保证
        能量和密度矩阵自洽；'printed' 用程序打印的最终能量。
    """

    # ---- 模型相关的默认值，子类 (Fe4S4) 覆盖 ----
    model_name = "Fe2S2"
    model_tag = "Fe2S2"
    norb_model = F2M_NORB
    dm1_name = F2M_DM1_NAME
    dm2_name = F2M_DM2_NAME
    default_nsegment = "0 0 10 10 0 0"
    default_mx_spintwo = 10
    exe_env = "F2M_ICIPT2_FE2S2"
    exe_name = "F2M_iCIPT2_Fe2S2.exe"

    def __init__(
        self,
        executable=None,
        workdir=None,
        spin=0,
        nroots=1,
        cmin=1e-4,
        etol=1e-10,
        mx_spintwo=None,
        mode=1,
        mx_ligand_ext_level=1,
        mx_ligand_spintwo=3,
        nsegment=None,
        nvalelec=2,
        irrep=0,
        selection=1,
        calc_rdm=True,
        dumprdm=None,
        energy_source="rdm",
        only_ferro_antiferro=False,
        hmat_brute_force_threshold=0,
        hmat_update=True,
        prune_init_guess=True,
        omp_num_threads=None,
        print_ci_each_call=True,
        n_top_csf=3,
        verbose=None,
    ):
        if mode not in (1, 2, 3):
            raise ValueError("mode must be 1, 2 or 3 (got %r)" % (mode,))
        if selection not in (0, 1):
            raise ValueError("selection must be 0 or 1 (got %r)" % (selection,))
        if executable is None:
            executable = _locate_executable(self.exe_env, self.exe_name, self.model_tag)
        self.executable = executable
        self.workdir = os.path.abspath(
            workdir if workdir is not None else "f2m_icipt2_work"
        )
        self.spin = spin
        self.nroots = nroots
        self.cmin = cmin
        self.etol = etol
        self.mx_spintwo = (
            self.default_mx_spintwo if mx_spintwo is None else int(mx_spintwo)
        )
        self.mode = int(mode)
        self.mx_ligand_ext_level = int(mx_ligand_ext_level)
        self.mx_ligand_spintwo = int(mx_ligand_spintwo)
        self.nsegment = self.default_nsegment if nsegment is None else nsegment
        self.nvalelec = nvalelec
        self.irrep = irrep
        self.selection = int(selection)
        self.calc_rdm = bool(calc_rdm)
        if dumprdm is None:
            dumprdm = 2 if self.calc_rdm else 0
        if int(dumprdm) not in (0, 1, 2, 7):
            raise ValueError("dumprdm must be 0, 1, 2 or 7 (got %r)" % (dumprdm,))
        self.dumprdm = int(dumprdm)
        if self.calc_rdm and self.dumprdm == 0:
            raise ValueError(
                "calc_rdm=True needs dumprdm = 1, 2 or 7 (dumprdm=0 writes no RDM)"
            )
        if energy_source not in ("rdm", "printed"):
            raise ValueError("energy_source must be 'rdm' or 'printed'")
        self.energy_source = energy_source
        self.only_ferro_antiferro = bool(only_ferro_antiferro)
        self.hmat_brute_force_threshold = int(hmat_brute_force_threshold)
        self.hmat_update = bool(hmat_update)
        self.prune_init_guess = bool(prune_init_guess)
        self.omp_num_threads = omp_num_threads
        self.print_ci_each_call = bool(print_ci_each_call)
        self.n_top_csf = int(n_top_csf)
        self.verbose = logger.NOTE if verbose is None else verbose

        self.converged = True
        self.e_ci = None
        self.e_printed = None
        self.mol = None

        # 诊断信息
        self.history = []  # 每次真正运行 F2M_iCIPT2 的记录
        self.last_output = None  # 最后一次运行的 stdout
        self.mc = None  # 由 driver 设置，供 print_summary 使用
        self.mc_history = []  # 宏迭代历史

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
        log.info("******** F2M_iCIPT2 solver (%s) ********", self.model_name)
        log.info("executable      = %s", self.executable)
        log.info("workdir         = %s", self.workdir)
        log.info("spin (2S)       = %d", self.spin)
        log.info("nroots          = %d", self.nroots)
        log.info("cmin / etol     = %s / %s", self.cmin, self.etol)
        log.info("mx_spintwo      = %d", self.mx_spintwo)
        log.info("mode            = %d", self.mode)
        log.info("mx_ligand_ext   = %d", self.mx_ligand_ext_level)
        log.info("mx_ligand_spin  = %d", self.mx_ligand_spintwo)
        log.info("nsegment        = %s", self.nsegment)
        log.info("selection       = %d", self.selection)
        log.info("dumprdm         = %d", self.dumprdm)
        log.info("only Fe/AF      = %d", int(self.only_ferro_antiferro))
        log.info("hmat brute thr  = %d", self.hmat_brute_force_threshold)
        log.info("hmat update     = %d", int(self.hmat_update))
        log.info("prune init guess= %d", int(self.prune_init_guess))
        log.info("energy source   = %s", self.energy_source)
        log.info(
            "OMP_NUM_THREADS = %s",
            self.omp_num_threads if self.omp_num_threads is not None else "(inherited)",
        )
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
            selection=self.selection,
            dumprdm=self.dumprdm,
            task=self._task_line(spintwo),
            spin=spintwo / 2.0,
            nstates=self.nroots,
        )
        with open(filename, "w") as f:
            f.write(text)
        return filename

    def _build_cmd(self, inputfile, integralfile):
        """F2M_iCIPT2 的命令行参数 (当前版本接口)。"""
        return [
            self.executable,
            inputfile,
            integralfile,
            str(self.mx_spintwo),
            str(self.mode),
            str(self.mx_ligand_ext_level),
            str(self.mx_ligand_spintwo),
            "1" if self.calc_rdm else "0",  # calculate_dm (程序已忽略)
            "1" if self.only_ferro_antiferro else "0",
            str(self.hmat_brute_force_threshold),
            "1" if self.hmat_update else "0",
            "1" if self.prune_init_guess else "0",
        ]

    def _run_program(self, inputfile, integralfile, spintwo):
        cmd = self._build_cmd(inputfile, integralfile)

        env = os.environ.copy()
        if self.omp_num_threads is not None:
            env["OMP_NUM_THREADS"] = str(int(self.omp_num_threads))

        logfile = os.path.join(
            self.workdir, "f2m_icipt2_%s.out" % self.model_name.lower()
        )
        errfile = os.path.join(
            self.workdir, "f2m_icipt2_%s.err" % self.model_name.lower()
        )
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
        if norb != self.norb_model:
            raise ValueError(
                "%s is compiled for the %d-orbital model, but norb = %d "
                "was passed.  All %d fcidump orbitals must be active."
                % (self.__class__.__name__, self.norb_model, norb, self.norb_model)
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
        for fname in (self.dm1_name, self.dm2_name):
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

        # 新版本程序在 "task N final energies:" 里打印最终能量；
        # 旧版本只有 Eig[...] 行，作为后备。
        energies = parse_task_energies(output, self.nroots)
        if energies is None:
            energies = parse_ci_energies(output, self.nroots)
        # FCIDUMP 里已经带了 ecore，所以打印出来的能量就是总能量
        self.e_printed = energies[0]

        if not self.calc_rdm:
            dm1 = dm2 = None
            e_rdm = None
        else:
            dm1_path = os.path.join(self.workdir, self.dm1_name)
            dm2_path = os.path.join(self.workdir, self.dm2_name)
            if not (os.path.isfile(dm1_path) and os.path.isfile(dm2_path)):
                raise RuntimeError(
                    "F2M_iCIPT2 did not dump the RDMs in %s (dumprdm=%d)"
                    % (self.workdir, self.dumprdm)
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
                        "Inconsistent 1-/2-RDM from F2M_iCIPT2 (err = %.3e)" % err
                    )

            e_rdm = rdm_energy(h1e, eri, dm1, dm2, ecore)

            if abs(e_rdm - self.e_printed) > 1e-6:
                logger.warn(
                    self,
                    "F2M_iCIPT2: printed energy (%.12f) differs from the RDM "
                    "energy (%.12f) by %.3e.  The RDM energy will be used.",
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

        # 诊断信息
        ncsf_list = _NCSF_RE.findall(output)
        ncsf = int(ncsf_list[-1]) if ncsf_list else None
        csfs = parse_important_csfs(output)
        self.last_output = output

        if self.print_ci_each_call and csfs:
            print(
                "\n[F2M_iCIPT2 %s call %d] E = %.12f, ncsf = %s, top CSFs:"
                % (self.model_name, self._ncall, e_ci, str(ncsf)),
                flush=True,
            )
            for c in csfs[: self.n_top_csf]:
                print(
                    "    %s   % .6f (%6.2f%%)" % (c["label"], c["coef"], c["weight"]),
                    flush=True,
                )

        self.history.append(
            dict(
                call=self._ncall,
                cmin=self.cmin,
                mode=self.mode,
                mx_spintwo=self.mx_spintwo,
                mx_ligand_ext_level=self.mx_ligand_ext_level,
                mx_ligand_spintwo=self.mx_ligand_spintwo,
                selection=self.selection,
                nelec=nelec_tot,
                spin_2S=spintwo,
                e_printed=self.e_printed,
                e_rdm=e_rdm,
                e_used=e_ci,
                ncsf=ncsf,
                csfs=csfs,
                logfile=logfile,
            )
        )
        return e_ci, ci

    def make_rdm1(self, fcivec=None, norb=None, nelec=None, **kwargs):
        if self._rdm1 is None:
            raise RuntimeError(
                "No RDM available.  Set calc_rdm=True (dumprdm=1/2/7) and call "
                "kernel() first."
            )
        return self._rdm1

    def make_rdm2(self, fcivec=None, norb=None, nelec=None, **kwargs):
        if self._rdm2 is None:
            raise RuntimeError(
                "No RDM available.  Set calc_rdm=True (dumprdm=2/7) and call "
                "kernel() first."
            )
        return self._rdm2

    def make_rdm12(self, fcivec=None, norb=None, nelec=None, **kwargs):
        return self.make_rdm1(), self.make_rdm2()

    def spin_square(self, fcivec=None, norb=None, nelec=None):
        s = float(self.spin) / 2.0
        return s * (s + 1.0), self.spin + 1

    #######################################################
    # 诊断输出
    #######################################################

    def dump_mo(self, mc=None, directory=None, basename="mo_coeff"):
        """把 mo_coeff 写成 csv，并画成热图。

        Returns
        -------
        (csv_path, png_path)
        """
        if mc is None:
            mc = self.mc
        if mc is None:
            raise RuntimeError("No MCSCF object available; call kernel() first.")
        if directory is None:
            directory = self.workdir
        os.makedirs(directory, exist_ok=True)

        mo = numpy.asarray(mc.mo_coeff)
        csv_path = os.path.join(directory, basename + ".csv")
        header = (
            "F2M_iCIPT2 %s: mo_coeff (rows: fcidump orbital, columns: MO); "
            "ncas=%d nelecas=%s E_tot=%.12f E_cas=%.12f"
            % (
                self.model_name,
                mc.ncas,
                str(mc.nelecas),
                mc.e_tot,
                mc.e_cas,
            )
        )
        numpy.savetxt(csv_path, mo, delimiter=",", fmt="%.10f", header=header)

        png_path = os.path.join(directory, basename + "_heatmap.png")
        try:
            plot_matrix_heatmap(
                mo,
                png_path,
                title="%s mo_coeff  (E = %.8f)" % (self.model_name, mc.e_tot),
                xlabel="MO index",
                ylabel="fcidump orbital index",
            )
        except ImportError:
            logger.warn(self, "matplotlib is not available; skip the mo_coeff heatmap")
            png_path = None
        return csv_path, png_path

    def print_mo(self, mc=None, stream=None):
        """打印 mo_coeff 矩阵。"""
        import sys

        if stream is None:
            stream = sys.stdout
        if mc is None:
            mc = self.mc
        if mc is None:
            raise RuntimeError("No MCSCF object available; call kernel() first.")
        stream.write(
            numpy.array2string(
                mc.mo_coeff, precision=6, suppress_small=True, max_line_width=250
            )
            + "\n"
        )
        stream.flush()

    def print_summary(self, mc=None, stream=None):
        """打印 MCSCF 迭代过程、最后一次 CI 调用的局域态/重要组态和能量。"""
        import sys

        if stream is None:
            stream = sys.stdout
        if mc is None:
            mc = self.mc

        def w(msg):
            stream.write(msg + "\n")

        line = "=" * 78
        w(line)
        w(" F2M_iCIPT2 (%s) summary" % self.model_name)
        w(line)
        w(" executable      : %s" % self.executable)
        w(" workdir         : %s" % self.workdir)
        w(" norb            : %d" % self.norb_model)
        w(" mode            : %d" % self.mode)
        w(" mx_spintwo      : %d" % self.mx_spintwo)
        w(
            " ligand ext/spin : %d / %d"
            % (self.mx_ligand_ext_level, self.mx_ligand_spintwo)
        )
        w(" selection       : %d" % self.selection)
        w(" dumprdm         : %d" % self.dumprdm)
        w(" only Fe/AF      : %d" % int(self.only_ferro_antiferro))
        w(
            " hmat brute thr  : %d, hmat update : %d, prune init : %d"
            % (
                self.hmat_brute_force_threshold,
                int(self.hmat_update),
                int(self.prune_init_guess),
            )
        )
        w(" energy source   : %s" % self.energy_source)
        w(" binary calls    : %d" % len(self.history))
        if mc is not None:
            nvirt = mc.mo_coeff.shape[1] - mc.ncore - mc.ncas
            w(" MCSCF converged : %s" % mc.converged)
            w(" MCSCF e_tot     : %.12f" % mc.e_tot)
            w(" MCSCF e_cas     : %.12f" % mc.e_cas)
            w(" ncore/ncas/nvirt: %d/%d/%d" % (mc.ncore, mc.ncas, nvirt))
            w(" mo_coeff:")
            for ln in numpy.array2string(
                mc.mo_coeff, precision=6, suppress_small=True, max_line_width=250
            ).splitlines():
                w("   " + ln)
            if getattr(mc, "mo_energy", None) is not None:
                w(" mo_energy:")
                w(
                    "   "
                    + numpy.array2string(
                        numpy.asarray(mc.mo_energy),
                        precision=6,
                        suppress_small=True,
                        max_line_width=250,
                    )
                )
        w("")

        if self.history:
            w(" CI solver calls (binary runs):")
            w(
                "   %-3s %-21s %-21s %-7s %s"
                % ("#", "E(printed)", "E(rdm)", "ncsf", "cmin")
            )
            for h in self.history:
                e_rdm = "%.12f" % h["e_rdm"] if h["e_rdm"] is not None else "N/A"
                w(
                    "   %-3d %-21.12f %-21s %-7s %s"
                    % (
                        h["call"],
                        h["e_printed"],
                        e_rdm,
                        str(h["ncsf"]),
                        str(h["cmin"]),
                    )
                )
            w("")

        if self.mc_history:
            w(" MCSCF macro iterations:")
            w(
                "   %-5s %-21s %-12s %-12s %-12s"
                % ("macro", "e_tot", "dE", "|g_orb|", "|ddm|")
            )
            for h in self.mc_history:
                w(
                    "   %-5s %-21s %-12s %-12s %-12s"
                    % (
                        str(h.get("imacro")),
                        "%.12f" % h["e_tot"] if h.get("e_tot") is not None else "N/A",
                        "%.3e" % h["de"] if h.get("de") is not None else "N/A",
                        (
                            "%.3e" % h["norm_gorb0"]
                            if h.get("norm_gorb0") is not None
                            else "N/A"
                        ),
                        (
                            "%.3e" % h["norm_ddm"]
                            if h.get("norm_ddm") is not None
                            else "N/A"
                        ),
                    )
                )
            w("")

        blocks = parse_local_states(self.last_output or "")
        if blocks:
            w(" Local states of the last F2M_iCIPT2 call:")
            for b in blocks:
                states = " ".join(
                    "(n=%d,2S=%d)x%d" % (st[0], st[1], st[3]) for st in b["states"]
                )
                w("   block %2d (norb=%s): %s" % (b["block"], str(b["norb"]), states))
            w("")

        csf = extract_important_csf(self.last_output or "")
        if csf:
            w(" Important CSFs of the last F2M_iCIPT2 call:")
            for ln in csf.splitlines():
                w("   " + ln)
            w("")

        for fname in ("mo_coeff.csv", "mo_coeff_heatmap.png"):
            fpath = os.path.join(self.workdir, fname)
            if os.path.isfile(fpath):
                w(" mo dump        : %s" % fpath)

        if self.history:
            last = self.history[-1]
            w(
                " Last call energy: E(printed) = %.12f  E(rdm) = %s  E(used) = %.12f"
                % (
                    last["e_printed"],
                    "%.12f" % last["e_rdm"] if last["e_rdm"] is not None else "N/A",
                    last["e_used"],
                )
            )
        w(line)
        stream.flush()


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
    mx_spintwo=None,
    mode=1,
    mx_ligand_ext_level=1,
    mx_ligand_spintwo=3,
    nsegment=None,
    nvalelec=2,
    selection=1,
    calc_rdm=True,
    dumprdm=None,
    energy_source="rdm",
    only_ferro_antiferro=False,
    hmat_brute_force_threshold=0,
    hmat_update=True,
    prune_init_guess=True,
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
    """MCSCF driver: FCIDUMP -> PySCF CASSCF，CI solver 用 F2M_iCIPT2。

    Parameters
    ----------
    fcidump : str
        模型 Hamiltonian 的 FCIDUMP 文件（例如
        ``Active-space-model-for-Iron-Sulfur-Clusters/Fe2S2_and_Fe4S4/Fe2S2/fe2s2``）。
    ncas, nelecas : int, optional
        活性轨道数和活性电子数。默认取 FCIDUMP 里的全部轨道 / 电子。
        F2M_iCIPT2 只能处理整个模型，所以 ``ncas`` 必须等于 FCIDUMP 的 NORB。
    active_orbs : None
        为了保持 F2M_iCIPT2 写死的轨道顺序，这里不允许重排轨道。
    spin : int
        2S，默认 0。
    cmin, etol, mx_spintwo, nsegment, nvalelec :
        透传给 F2M_iCIPT2（见 :class:`F2MiCIPT2Solver`）。
    mode : int
        1: 只有自旋激发 (Fe2+/Fe3+)；2: + MMCT (Fe+/Fe4+)；
        3: + LMCT (配体激发)。
    mx_ligand_ext_level, mx_ligand_spintwo :
        mode 3 下配体局域激发的等级和最大局域自旋。
    selection : int
        0 = 不做 selection；1 = 做 selection（由 cmin 控制）。
    calc_rdm, dumprdm :
        是否请求 RDM 以及 ``dumprdm`` 的取值（0/1/2/7）。
    only_ferro_antiferro, hmat_brute_force_threshold, hmat_update,
    prune_init_guess :
        透传给 F2M_iCIPT2 的新选项。
    energy_source : {'rdm', 'printed'}
        透传给 F2M_iCIPT2。
    omp_num_threads : int, optional
        给 F2M_iCIPT2 子进程设置的 OMP_NUM_THREADS；None 表示继承环境。
    run_mcscf : bool
        True  : 走完整的 CASSCF 宏迭代；
        False : 只做一次 CASCI（一次 F2M_iCIPT2 调用）。
    internal_rotation : bool
        是否允许活性空间内部的轨道旋转。
    print_summary : bool
        结束时打印宏迭代历史、最后一次 F2M_iCIPT2 调用的局域态和重要组态、
        能量以及 mo_coeff。
    dump_mo : bool
        结束时把 mo_coeff 写成 csv 并画热图。

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
            "F2M_iCIPT2 requires all %d fcidump orbitals to be active "
            "(ncas = %d was requested).  Its local-orbital structure and "
            "ref_occ are hard-coded for the full model." % (norb, ncas)
        )
    if active_orbs is not None:
        raise ValueError(
            "active_orbs is not supported: F2M_iCIPT2 requires the original "
            "fcidump orbital order."
        )
    norb_model = solver.norb_model if solver is not None else F2M_NORB
    if norb != norb_model:
        raise ValueError(
            "F2M_iCIPT2 for this model is compiled for %d orbitals, "
            "this fcidump has %d" % (norb_model, norb)
        )
    if run_mcscf and not calc_rdm:
        raise ValueError("MCSCF needs the RDMs of the CI solver; set calc_rdm=True.")

    log = logger.new_logger(None, verbose)
    log.info("")
    log.info("******** MCSCF with F2M_iCIPT2 solver ********")
    log.info("fcidump = %s", fcidump)
    log.info(
        "norb    = %d, nel = %d, ncas = %d, nelecas = %s",
        norb,
        nelec_fcidump,
        ncas,
        str(nelecas),
    )

    mol, mf, mc = setup_mcscf_from_fcidump(
        fcidump,
        ncas=ncas,
        nelecas=nelecas,
        active_orbs=None,
        mc_type="casscf",
        verbose=verbose,
    )

    # 关键: 初猜使用 FCIDUMP 的轨道基。
    mo_init = numpy.eye(norb)
    mc.mo_coeff = mo_init

    if solver is None:
        solver = F2MiCIPT2Solver(
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
            selection=selection,
            calc_rdm=calc_rdm,
            dumprdm=dumprdm,
            energy_source=energy_source,
            only_ferro_antiferro=only_ferro_antiferro,
            hmat_brute_force_threshold=hmat_brute_force_threshold,
            hmat_update=hmat_update,
            prune_init_guess=prune_init_guess,
            omp_num_threads=omp_num_threads,
            print_ci_each_call=print_ci_each_call,
            n_top_csf=n_top_csf,
            verbose=verbose,
        )
    mc.fcisolver = solver
    solver.mc = mc

    mc.conv_tol = mc_conv_tol
    mc.max_cycle_macro = mc_max_macro
    mc.internal_rotation = internal_rotation
    mc.canonicalization = False

    macro_history = []

    def _record_macro(envs):
        try:
            macro_history.append(
                dict(
                    imacro=envs.get("imacro"),
                    e_tot=envs.get("e_tot"),
                    e_cas=envs.get("e_cas"),
                    de=envs.get("de"),
                    norm_gorb0=envs.get("norm_gorb0"),
                    norm_ddm=envs.get("norm_ddm"),
                )
            )
        except Exception:
            pass

    if run_mcscf:
        mc.callback = _record_macro
        mc.kernel(mo_init)

        mo_final = mc.mo_coeff
        offdiag = numpy.abs(mo_final) - numpy.abs(numpy.eye(norb))
        if offdiag.max() > 1e-6:
            logger.info(
                mc,
                "mo_coeff was rotated during MCSCF (max |offdiag| = %.3e). "
                "The F2M_iCIPT2 local-state model is defined in the FCIDUMP "
                "orbital basis, so the result depends on the basis.",
                offdiag.max(),
            )
    else:
        mc.e_tot, mc.e_cas, mc.ci = mc.casci(mo_init)

    if macro_history:
        dedup = {}
        for h in macro_history:
            if h.get("imacro") is not None:
                dedup[h["imacro"]] = h
        solver.mc_history = [dedup[k] for k in sorted(dedup)]

    log.info("")
    log.info("MCSCF/F2M_iCIPT2 total energy = %.12f", mc.e_tot)
    log.info("MCSCF/F2M_iCIPT2 CI    energy = %.12f", mc.e_cas)

    if dump_mo:
        csv_path, png_path = solver.dump_mo(mc)
        log.info("mo_coeff dumped to %s", csv_path)
        if png_path is not None:
            log.info("mo_coeff heatmap: %s", png_path)

    if print_summary:
        solver.print_summary(mc=mc)

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

    mc = kernel(fcidump, cmin=1e-1, etol=1e-4, mode=1, verbose=4, print_summary=True)
    print("E = ", mc.e_tot)
