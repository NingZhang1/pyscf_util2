"""解析 cNEO-iCIPT2 输出文件（ES / PES 两套结果）。

典型输出格式::

    *******************************************************
    *       cNEO-iCIPT2 ES       *
    -------------------------------------------------------
    Space : spintwo 0 irrep 0
    ------------------------------------------------------------------------------
    |       cmin |       ncsf |       ncfg | (   eiCI,   ePT2,   etot) |
    ------------------------------------------------------------------------------
    | 1.0000e-04 |      10908 |       7823 | (  -99.247548998579,  ...) |
    ...
    -------------------------------------------------------
    Nuclear Kinetic Energy
    State 0
    1.0000e-04 |       0.468578994258
    ...
    -------------------------------------------------------
    *       cNEO-iCIPT2 PES      *
    -------------------------------------------------------
    Space : spintwo 0 irrep 0
    ...

一个文件里可能出现多组 ``ES`` / ``PES`` 块；每个块内可以有多个
``Space : spintwo ... irrep ...`` 扇区；每个 cmin 行可以带一个或多个
``(eiCI, ePT2, etot)`` 根（每个根对应一个电子态）。ES 块还额外带有
``Nuclear Kinetic Energy`` 段，按 ``State N`` 给出每个态的核动能。

本模块把上述内容解析为结构化对象，并提供“按 ES / PES 分别统计”的接口。

用法::

    from pyscf_util.Analyzer.iCIPT2.analyzer_cneo import (
        read_cneo_icipt2, extract_cneo_states, print_cneo_summary,
    )

    blocks = read_cneo_icipt2("cneo.out")
    print_cneo_summary(blocks)              # 先 ES 后 PES，各自成表

    es_states = extract_cneo_states("cneo.out", kind="ES")
    pes_states = extract_cneo_states("cneo.out", kind="PES")

命令行::

    python -m pyscf_util.Analyzer.iCIPT2.analyzer_cneo cneo.out
    python -m pyscf_util.Analyzer.iCIPT2.analyzer_cneo cneo.out --kind ES
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass, field
from typing import Iterable, Iterator, Optional

from tabulate import tabulate

__all__ = [
    "CneoRoot",
    "CneoRow",
    "CneoSpace",
    "CneoBlock",
    "CneoState",
    "parse_cneo_icipt2",
    "read_cneo_icipt2",
    "iter_cneo_states",
    "extract_cneo_states",
    "split_by_kind",
    "cneo_statistics",
    "print_cneo_states",
    "print_cneo_summary",
]

#: 块标题里的关键字
TITLE = "cNEO-iCIPT2"

_FLOAT = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"
_FLOAT_RE = re.compile(rf"^{_FLOAT}$")
_ROOT_RE = re.compile(rf"\(\s*({_FLOAT})\s*,\s*({_FLOAT})\s*,\s*({_FLOAT})\s*\)")
_HEADER_RE = re.compile(r"cNEO-iCIPT2\s+(ES|PES)\b", re.IGNORECASE)
_SPACE_RE = re.compile(r"spintwo\s+(-?\d+)\s+irrep\s+(-?\d+)")
_STATE_RE = re.compile(r"^State\s+(\d+)")
_END_RE = re.compile(r"^\s*\*{4,}\s*$")


# ---------------------------------------------------------------------------
# 数据结构
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CneoRoot:
    """一个 ``(eiCI, ePT2, etot)`` 能量三元组，对应一个电子态。"""

    eiCI: float
    ePT2: float
    etot: float


@dataclass(frozen=True)
class CneoRow:
    """表格中的一行：一个 cmin 及其在该 cmin 下的一个或多个根。"""

    cmin: float
    ncsf: int
    ncfg: int
    roots: tuple[CneoRoot, ...]


@dataclass
class CneoSpace:
    """一个 ``Space : spintwo X irrep Y`` 扇区。"""

    block: int
    kind: str
    spintwo: int
    irrep: int
    rows: list[CneoRow] = field(default_factory=list)
    #: 核动能：state -> [(cmin, T_nuc), ...]，仅 ES 块有
    nuclear_kinetic: dict[int, list[tuple[float, float]]] = field(
        default_factory=dict
    )
    line_number: int = 0

    @property
    def nstates(self) -> int:
        """该扇区的态数（以每个 cmin 行最多的根数为准）。"""
        return max((len(row.roots) for row in self.rows), default=0)

    def nuclear_kinetic_map(self, state: int) -> dict[float, float]:
        """返回某个态的 ``cmin -> T_nuc`` 映射。"""
        return dict(self.nuclear_kinetic.get(state, []))


@dataclass
class CneoBlock:
    """一个 ``cNEO-iCIPT2 ES`` 或 ``cNEO-iCIPT2 PES`` 块。"""

    index: int
    kind: str
    spaces: list[CneoSpace] = field(default_factory=list)
    line_number: int = 0

    def __len__(self) -> int:
        return len(self.spaces)

    def __iter__(self) -> Iterator[CneoSpace]:
        return iter(self.spaces)

    def __getitem__(self, item):
        return self.spaces[item]


@dataclass(frozen=True)
class CneoState:
    """扁平化后的“一个态在一个 cmin 下的结果”。

    这是“按态统计”的最基本单位：``(kind, block, spintwo, irrep, state, cmin)``。
    """

    kind: str
    block: int
    spintwo: int
    irrep: int
    state: int
    cmin: float
    ncsf: int
    ncfg: int
    eiCI: float
    ePT2: float
    etot: float
    nuclear_kinetic: Optional[float] = None
    line_number: int = 0

    def as_dict(self) -> dict:
        return {
            "kind": self.kind,
            "block": self.block,
            "spintwo": self.spintwo,
            "irrep": self.irrep,
            "state": self.state,
            "cmin": self.cmin,
            "ncsf": self.ncsf,
            "ncfg": self.ncfg,
            "eiCI": self.eiCI,
            "ePT2": self.ePT2,
            "etot": self.etot,
            "nuclear_kinetic": self.nuclear_kinetic,
            "line_number": self.line_number,
        }

    def __str__(self) -> str:
        nk = "" if self.nuclear_kinetic is None else f", T_nuc {self.nuclear_kinetic:.12f}"
        return (
            f"{self.kind} block {self.block}, spintwo {self.spintwo}, "
            f"irrep {self.irrep}, state {self.state}, cmin {self.cmin:.4e}: "
            f"ncsf {self.ncsf}, ncfg {self.ncfg}, "
            f"eiCI {self.eiCI:.12f}, ePT2 {self.ePT2:.12f}, "
            f"etot {self.etot:.12f}{nk}"
        )


# ---------------------------------------------------------------------------
# 小工具
# ---------------------------------------------------------------------------


def _is_float(text: str) -> bool:
    return bool(_FLOAT_RE.match(text))


def _to_float(text: str) -> float:
    # 兼容 Fortran 的 D 指数（如 1.0D-03）
    return float(text.replace("D", "E").replace("d", "e"))


def _to_int(text: str) -> Optional[int]:
    try:
        return int(text)
    except (TypeError, ValueError):
        return None


def _split_pipe(line: str) -> list[str]:
    stripped = line.strip()
    if stripped.startswith("|"):
        stripped = stripped[1:]
    if stripped.endswith("|"):
        stripped = stripped[:-1]
    return [part.strip() for part in stripped.split("|")]


def _parse_table_row(line: str) -> Optional[CneoRow]:
    """解析 ``| cmin | ncsf | ncfg | (eiCI, ePT2, etot) ... |`` 一行。"""
    parts = _split_pipe(line)
    if len(parts) != 4:
        return None
    if not _is_float(parts[0]):
        return None  # 表头行
    cmin = _to_float(parts[0])
    ncsf = _to_int(parts[1])
    ncfg = _to_int(parts[2])
    if ncsf is None or ncfg is None:
        return None
    roots = tuple(
        CneoRoot(*(_to_float(g) for g in match.groups()))
        for match in _ROOT_RE.finditer(parts[3])
    )
    if not roots:
        return None
    return CneoRow(cmin=cmin, ncsf=ncsf, ncfg=ncfg, roots=roots)


def _parse_nuclear_kinetic_row(line: str) -> Optional[tuple[float, float]]:
    """解析核动能段的一行 ``cmin | T_nuc``。"""
    parts = _split_pipe(line)
    if len(parts) != 2 or not all(_is_float(p) for p in parts):
        return None
    return _to_float(parts[0]), _to_float(parts[1])


# ---------------------------------------------------------------------------
# 解析
# ---------------------------------------------------------------------------


def parse_cneo_icipt2(text: str) -> list[CneoBlock]:
    """解析文本中所有 cNEO-iCIPT2 的 ES / PES 块。

    Args:
        text: 文件内容。

    Returns:
        ``CneoBlock`` 列表，按文件出现顺序；每个块内 ``spaces`` 按出现顺序。
        没有找到任何块时返回空列表。
    """
    lines = text.splitlines()
    blocks: list[CneoBlock] = []

    block: Optional[CneoBlock] = None
    space: Optional[CneoSpace] = None
    in_nuclear_kinetic = False
    nuclear_state: Optional[int] = None

    def close_space() -> None:
        nonlocal space
        if space is not None and block is not None:
            block.spaces.append(space)
        space = None

    def close_block() -> None:
        nonlocal block
        close_space()
        if block is not None:
            blocks.append(block)
        block = None

    for line_number, line in enumerate(lines, start=1):
        header = _HEADER_RE.search(line)
        if header is not None:
            close_block()
            block = CneoBlock(
                index=len(blocks),
                kind=header.group(1).upper(),
                line_number=line_number,
            )
            in_nuclear_kinetic = False
            nuclear_state = None
            continue

        if _END_RE.match(line):
            close_block()
            in_nuclear_kinetic = False
            nuclear_state = None
            continue

        if block is None:
            continue

        match_space = _SPACE_RE.search(line)
        if match_space is not None:
            close_space()
            space = CneoSpace(
                block=block.index,
                kind=block.kind,
                spintwo=int(match_space.group(1)),
                irrep=int(match_space.group(2)),
                line_number=line_number,
            )
            in_nuclear_kinetic = False
            nuclear_state = None
            continue

        if space is None:
            continue

        if "Nuclear Kinetic Energy" in line:
            in_nuclear_kinetic = True
            nuclear_state = None
            continue

        is_table_row = line.lstrip().startswith("|")

        if in_nuclear_kinetic and not is_table_row:
            if line.strip().startswith("-"):
                in_nuclear_kinetic = False
                nuclear_state = None
                continue
            match_state = _STATE_RE.match(line.strip())
            if match_state is not None:
                nuclear_state = int(match_state.group(1))
                continue
            pair = _parse_nuclear_kinetic_row(line)
            if pair is not None:
                space.nuclear_kinetic.setdefault(
                    0 if nuclear_state is None else nuclear_state, []
                ).append(pair)
            continue

        if in_nuclear_kinetic:
            # 缺少分隔线时，遇到表格行即认为核动能段结束
            in_nuclear_kinetic = False
            nuclear_state = None

        if is_table_row:
            row = _parse_table_row(line)
            if row is not None:
                space.rows.append(row)

    close_block()
    return blocks


def read_cneo_icipt2(
    filename: str, encoding: str = "utf-8"
) -> list[CneoBlock]:
    """读取文件并解析其中所有 cNEO-iCIPT2 块。"""
    with open(filename, "r", encoding=encoding, errors="replace") as handle:
        return parse_cneo_icipt2(handle.read())


# ---------------------------------------------------------------------------
# 按态展开 / 按 ES、PES 分组
# ---------------------------------------------------------------------------


def iter_cneo_states(blocks: Iterable[CneoBlock]) -> Iterator[CneoState]:
    """把块展开为“一行一个态一个 cmin”的 ``CneoState`` 序列。

    同一 cmin 行若有多个根，按根的顺序编号为 ``state 0, 1, ...``；ES 块中
    对应 ``State N`` 的核动能会填到 ``nuclear_kinetic``。
    """
    for block in blocks:
        for space in block.spaces:
            nuclear_maps = {
                state: space.nuclear_kinetic_map(state)
                for state in space.nuclear_kinetic
            }
            for row in space.rows:
                for state, root in enumerate(row.roots):
                    yield CneoState(
                        kind=block.kind,
                        block=block.index,
                        spintwo=space.spintwo,
                        irrep=space.irrep,
                        state=state,
                        cmin=row.cmin,
                        ncsf=row.ncsf,
                        ncfg=row.ncfg,
                        eiCI=root.eiCI,
                        ePT2=root.ePT2,
                        etot=root.etot,
                        nuclear_kinetic=nuclear_maps.get(state, {}).get(row.cmin),
                        line_number=space.line_number,
                    )


def extract_cneo_states(
    filename: str,
    kind: Optional[str] = None,
    encoding: str = "utf-8",
) -> list[CneoState]:
    """读取文件并返回按态展开的结果。

    Args:
        filename: 文件路径。
        kind: ``"ES"``、``"PES"`` 或 ``None``（两者都要）。
    """
    blocks = read_cneo_icipt2(filename, encoding)
    if kind is not None:
        kind = kind.upper()
        blocks = [block for block in blocks if block.kind == kind]
    return list(iter_cneo_states(blocks))


def split_by_kind(blocks: Iterable[CneoBlock]) -> dict[str, list[CneoBlock]]:
    """按 ``ES`` / ``PES`` 分组，返回有序 dict。"""
    grouped: dict[str, list[CneoBlock]] = {}
    for block in blocks:
        grouped.setdefault(block.kind, []).append(block)
    return grouped


def cneo_statistics(blocks: Iterable[CneoBlock]) -> dict[str, dict[str, int]]:
    """统计每种块（ES/PES）的块数、扇区数、电子态数与记录数。

    返回的每个字典包含：

    - ``blocks``: 块数
    - ``spaces``: ``Space :`` 扇区数
    - ``states``: 电子态数（每个扇区取各 cmin 行根数的最大值后求和）
    - ``cmin_points``: cmin 点数（表格行数）
    - ``records``: 数据记录数（cmin 点数 × 每行根数）
    """
    stats: dict[str, dict[str, int]] = {}
    for block in blocks:
        entry = stats.setdefault(
            block.kind,
            {
                "blocks": 0,
                "spaces": 0,
                "states": 0,
                "cmin_points": 0,
                "records": 0,
            },
        )
        entry["blocks"] += 1
        entry["spaces"] += len(block.spaces)
        for space in block.spaces:
            entry["cmin_points"] += len(space.rows)
            entry["states"] += space.nstates
            entry["records"] += sum(len(row.roots) for row in space.rows)
    return stats


# ---------------------------------------------------------------------------
# 输出
# ---------------------------------------------------------------------------


def print_cneo_states(
    states: Iterable[CneoState],
    tablefmt: str = "grid",
    show_nuclear_kinetic: Optional[bool] = None,
) -> None:
    """把 ``CneoState`` 列表打印成表格（一行一个态）。"""
    states = list(states)
    if show_nuclear_kinetic is None:
        show_nuclear_kinetic = any(
            state.nuclear_kinetic is not None for state in states
        )

    headers = [
        "kind",
        "block",
        "spintwo",
        "irrep",
        "state",
        "cmin",
        "ncsf",
        "ncfg",
        "eiCI",
        "ePT2",
        "etot",
    ]
    rows = []
    for state in states:
        row = [
            state.kind,
            state.block,
            state.spintwo,
            state.irrep,
            state.state,
            f"{state.cmin:.4e}",
            state.ncsf,
            state.ncfg,
            state.eiCI,
            state.ePT2,
            state.etot,
        ]
        if show_nuclear_kinetic:
            row.append(state.nuclear_kinetic)
        rows.append(row)
    if show_nuclear_kinetic:
        headers.append("T_nuc")

    print(
        tabulate(
            rows,
            headers,
            tablefmt=tablefmt,
            stralign="right",
            numalign="right",
            floatfmt=".12f",
            missingval="-",
            # cmin 列保持 "1.0000e-04" 的科学计数法，不让 tabulate 重新解析
            disable_numparse=[5],
        )
    )


def print_cneo_summary(
    blocks: Iterable[CneoBlock],
    tablefmt: str = "grid",
) -> None:
    """分别打印 ES 与 PES 的结果（各自一张表），并给出统计行。"""
    blocks = list(blocks)
    grouped = split_by_kind(blocks)
    stats = cneo_statistics(blocks)

    for kind, kind_blocks in grouped.items():
        info = stats[kind]
        print(
            f"=== cNEO-iCIPT2 {kind}: "
            f"{info['blocks']} block(s), {info['spaces']} space(s), "
            f"{info['states']} state(s), {info['cmin_points']} cmin point(s), "
            f"{info['records']} record(s) ==="
        )
        print_cneo_states(iter_cneo_states(kind_blocks), tablefmt=tablefmt)
        print()


# ---------------------------------------------------------------------------
# 命令行入口
# ---------------------------------------------------------------------------


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="解析 cNEO-iCIPT2 输出，分别统计 ES / PES 结果。",
    )
    parser.add_argument("files", nargs="+", help="cNEO-iCIPT2 输出文件")
    parser.add_argument(
        "--kind",
        choices=["ES", "PES", "es", "pes", "all"],
        default="all",
        help="只输出 ES 或 PES（默认两者）",
    )
    parser.add_argument(
        "--tablefmt",
        default="grid",
        help="tabulate 表格格式（默认 grid，可选 plain/tsv/...）",
    )
    args = parser.parse_args(argv)

    kind = None if args.kind.lower() == "all" else args.kind.upper()

    for path in args.files:
        blocks = read_cneo_icipt2(path)
        if len(args.files) > 1:
            print(f"=== {path} ===")
        if not blocks:
            print(f"未在 {path} 中找到 cNEO-iCIPT2 块。")
            continue
        if kind is not None:
            blocks = [block for block in blocks if block.kind == kind]
            if not blocks:
                print(f"文件 {path} 中没有 {kind} 块。")
                continue
        print_cneo_summary(blocks, tablefmt=args.tablefmt)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
