"""解析 iCIPT2 输出文件里的 ``iCI_ENPT(2)_NonRela::Info`` 代码块。

iCIPT2 在不同版本/不同运行模式下会打印两种表头：

1. 多态格式（带自旋与对称性）::

    *************** iCI_ENPT(2)_NonRela::Info ***************

    ________________________________________________________________________
       spintwo |   symmetry |       ncfg |       ncsf |        evar |  ...
    ________________________________________________________________________
             0 |          0 |      95932 |     217314 | -75.45475170 |  ...
               |            |            |            | -75.26317999 |  ...
               |            |            |            | -75.20765095 |  ...
             0 |          1 |      52390 |     114358 | -75.20760436 |  ...
    ________________________________________________________________________

    *******************************************************

2. 旧格式（单态）::

    *************** iCI_ENPT(2)_NonRela::Info ***************

    __________________________________________________
          ncfg |       ncsf |        evar |   ept |  ...
    __________________________________________________
        129890 |     231392 | -547.69018534 |  ...
    __________________________________________________

    *******************************************************

本模块把两种格式统一解析为“按态（state）”存放的数据结构：每个数据行就是一个
电子态；同一 ``(spintwo, symmetry)`` 扇区内的后续态会把前四列留空，只写
``evar | ept | etot``，解析时自动继承上一行的取值。一个文件里可能出现多个
代码块（每个迭代步一个），``Enpt2State.block`` 用来区分它们。

用法::

    from pyscf_util.Analyzer.iCIPT2.enpt2_info import (
        extract_enpt2_states, group_states_by_sector, print_enpt2_states,
    )

    states = extract_enpt2_states("iCI.out")
    print_enpt2_states(states)
    for (spintwo, symmetry), sub in group_states_by_sector(states).items():
        print(spintwo, symmetry, [s.etot for s in sub])

命令行::

    python -m pyscf_util.Analyzer.iCIPT2.enpt2_info iCI.out
    python -m pyscf_util.Analyzer.iCIPT2.enpt2_info iCI.out --sector
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass, field
from typing import Iterable, Iterator, Optional

from tabulate import tabulate

__all__ = [
    "Enpt2State",
    "Enpt2Block",
    "parse_enpt2_info",
    "read_enpt2_info",
    "extract_enpt2_states",
    "group_states_by_sector",
    "print_enpt2_states",
]

#: 代码块标题里用于识别的关键字
BLOCK_TITLE = "iCI_ENPT(2)_NonRela::Info"

#: 表格分隔线，例如 "______"
_RULE_RE = re.compile(r"^\s*_{5,}\s*$")
#: 代码块结束的星号行，例如 "**********"
_END_RE = re.compile(r"^\s*\*{4,}\s*$")
#: 浮点数，兼容 Fortran 的 D 指数（如 1.0D-03）
_FLOAT_RE = re.compile(r"^[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?$")


# ---------------------------------------------------------------------------
# 数据结构
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Enpt2State:
    """``iCI_ENPT(2)_NonRela::Info`` 表格中的一行，即一个电子态。

    Attributes:
        block: 代码块序号（从 0 开始），对应文件里第几个 ``iCI_ENPT(2)`` 块。
        spintwo: 自旋（2S），旧格式没有该列时为 ``None``。
        symmetry: 空间对称性（irrep 编号），旧格式没有该列时为 ``None``。
        state: 该 ``(spintwo, symmetry)`` 扇区内第几个态，从 0 开始。
        ncfg: 组态数。
        ncsf: CSF 数。
        evar: 变分能量（iCI 能量）。
        ept: 二阶微扰修正。
        etot: 总能量（evar + ept）。
        line_number: 该行在源文件中的行号（从 1 开始）。
    """

    block: int
    spintwo: Optional[int]
    symmetry: Optional[int]
    state: int
    ncfg: int
    ncsf: int
    evar: float
    ept: float
    etot: float
    line_number: int = 0

    @property
    def sector(self) -> tuple[Optional[int], Optional[int]]:
        """返回 ``(spintwo, symmetry)``。"""
        return (self.spintwo, self.symmetry)

    def as_dict(self) -> dict:
        return {
            "block": self.block,
            "spintwo": self.spintwo,
            "symmetry": self.symmetry,
            "state": self.state,
            "ncfg": self.ncfg,
            "ncsf": self.ncsf,
            "evar": self.evar,
            "ept": self.ept,
            "etot": self.etot,
            "line_number": self.line_number,
        }

    def __str__(self) -> str:
        spin = "-" if self.spintwo is None else self.spintwo
        sym = "-" if self.symmetry is None else self.symmetry
        return (
            f"block {self.block}, spintwo {spin}, symmetry {sym}, state {self.state}: "
            f"ncfg {self.ncfg}, ncsf {self.ncsf}, "
            f"evar {self.evar:.8f}, ept {self.ept:.8f}, etot {self.etot:.8f}"
        )


@dataclass
class Enpt2Block:
    """文件中的一个 ``iCI_ENPT(2)_NonRela::Info`` 代码块。

    Attributes:
        index: 代码块序号（从 0 开始）。
        has_spin_symmetry: 表头是否带 ``spintwo | symmetry`` 两列。
        states: 该块内解析出的所有态（按文件顺序）。
        line_number: 块标题所在行号（从 1 开始）。
    """

    index: int
    has_spin_symmetry: bool
    states: list[Enpt2State] = field(default_factory=list)
    line_number: int = 0

    def __len__(self) -> int:
        return len(self.states)

    def __iter__(self) -> Iterator[Enpt2State]:
        return iter(self.states)

    def __getitem__(self, item):
        return self.states[item]

    def by_sector(self) -> dict[tuple[Optional[int], Optional[int]], list[Enpt2State]]:
        """按 ``(spintwo, symmetry)`` 分组的态。"""
        return group_states_by_sector(self.states)


# ---------------------------------------------------------------------------
# 解析
# ---------------------------------------------------------------------------


def _split_fields(line: str) -> list[str]:
    """按 ``|`` 切分一行，去掉两端空白和行尾的空字段。"""
    if "|" not in line:
        return []
    fields = [f.strip() for f in line.split("|")]
    while fields and fields[-1] == "":
        fields.pop()
    return fields


def _is_float(text: str) -> bool:
    return bool(_FLOAT_RE.match(text))


def _to_float(text: str) -> float:
    # Fortran 的 D 指数（如 1.0D-03）Python 不认识，替换成 E
    return float(text.replace("D", "E").replace("d", "e"))


def _to_int(text: str) -> Optional[int]:
    try:
        return int(text)
    except (TypeError, ValueError):
        return None


def _is_header(fields: Iterable[str]) -> bool:
    names = {f.lower() for f in fields}
    return {"evar", "ept", "etot"} <= names


def _parse_data_line(
    line: str,
    has_spin_symmetry: bool,
    last: dict,
) -> Optional[tuple]:
    """解析一行数据。

    Args:
        line: 原始文本行。
        has_spin_symmetry: 表头是否有 spintwo/symmetry 两列。
        last: 上一行解析出的字段，用于补全留空的列；成功时会被就地更新。

    Returns:
        ``(spintwo, symmetry, ncfg, ncsf, evar, ept, etot)``；不是数据行时返回
        ``None``。
    """
    fields = _split_fields(line)
    if len(fields) < 5:
        return None
    if not all(_is_float(f) for f in fields[-3:]):
        return None

    evar, ept, etot = (_to_float(f) for f in fields[-3:])
    lead = fields[:-3]

    new = dict(last)

    def _take_int(offset: int, key: str) -> None:
        value = _to_int(lead[offset])
        if value is not None:
            new[key] = value

    if has_spin_symmetry:
        if len(lead) == 4:
            _take_int(0, "spintwo")
            _take_int(1, "symmetry")
            _take_int(2, "ncfg")
            _take_int(3, "ncsf")
        elif len(lead) == 2:
            # 自旋/对称性两列被整体省略，继承上一行
            _take_int(0, "ncfg")
            _take_int(1, "ncsf")
        elif len(lead) == 0:
            pass  # 全部继承
        else:
            return None
    else:
        if len(lead) == 2:
            _take_int(0, "ncfg")
            _take_int(1, "ncsf")
        elif len(lead) == 0:
            pass
        else:
            return None

    if new["ncfg"] is None or new["ncsf"] is None:
        return None

    last.update(new)
    return (
        new["spintwo"],
        new["symmetry"],
        new["ncfg"],
        new["ncsf"],
        evar,
        ept,
        etot,
    )


def parse_enpt2_info(text: str) -> list[Enpt2Block]:
    """解析文本中所有 ``iCI_ENPT(2)_NonRela::Info`` 代码块。

    Args:
        text: 文件内容（字符串）。

    Returns:
        ``Enpt2Block`` 列表，按在文件中出现的顺序排列；每个块内的 ``states``
        按文件顺序排列。没有找到任何代码块时返回空列表。
    """
    lines = text.splitlines()
    blocks: list[Enpt2Block] = []
    index = 0
    block_index = -1
    total = len(lines)

    while index < total:
        if BLOCK_TITLE not in lines[index]:
            index += 1
            continue

        block_index += 1
        start_line = index + 1
        index += 1

        # --- 定位表头 ---
        header_fields: Optional[list[str]] = None
        while index < total:
            line = lines[index]
            if BLOCK_TITLE in line or _END_RE.match(line):
                break
            fields = _split_fields(line)
            if fields and _is_header(fields):
                header_fields = fields
                index += 1
                break
            index += 1

        if header_fields is None:
            # 只有标题、没有表头，跳过
            continue

        has_spin_symmetry = any(f.lower() == "spintwo" for f in header_fields)
        block = Enpt2Block(
            index=block_index,
            has_spin_symmetry=has_spin_symmetry,
            line_number=start_line,
        )

        # --- 解析数据行 ---
        last = {"spintwo": None, "symmetry": None, "ncfg": None, "ncsf": None}
        sector_counter: dict[tuple, int] = {}

        while index < total:
            line = lines[index]
            if BLOCK_TITLE in line or _END_RE.match(line):
                break
            if _RULE_RE.match(line) or not line.strip():
                index += 1
                continue

            parsed = _parse_data_line(line, has_spin_symmetry, last)
            if parsed is not None:
                spintwo, symmetry, ncfg, ncsf, evar, ept, etot = parsed
                sector = (spintwo, symmetry)
                state = sector_counter.get(sector, 0)
                sector_counter[sector] = state + 1
                block.states.append(
                    Enpt2State(
                        block=block_index,
                        spintwo=spintwo,
                        symmetry=symmetry,
                        state=state,
                        ncfg=ncfg,
                        ncsf=ncsf,
                        evar=evar,
                        ept=ept,
                        etot=etot,
                        line_number=index + 1,
                    )
                )
            index += 1

        blocks.append(block)

    return blocks


def read_enpt2_info(filename: str, encoding: str = "utf-8") -> list[Enpt2Block]:
    """读取文件并解析其中所有 ``iCI_ENPT(2)_NonRela::Info`` 代码块。"""
    with open(filename, "r", encoding=encoding, errors="replace") as handle:
        return parse_enpt2_info(handle.read())


def extract_enpt2_states(
    filename: str, encoding: str = "utf-8"
) -> list[Enpt2State]:
    """读取文件，返回所有代码块中的所有态（按文件顺序，扁平列表）。"""
    return [state for block in read_enpt2_info(filename, encoding) for state in block]


# ---------------------------------------------------------------------------
# 输出
# ---------------------------------------------------------------------------


def group_states_by_sector(
    states: Iterable[Enpt2State],
    per_block: bool = False,
) -> dict[tuple, list[Enpt2State]]:
    """按 ``(spintwo, symmetry)`` 对态分组。

    Args:
        states: ``Enpt2State`` 可迭代对象。
        per_block: 为 ``True`` 时键为 ``(block, spintwo, symmetry)``，即同一
            扇区在不同迭代步的态分开存放；默认合并。

    Returns:
        有序字典（普通 dict，保持首次出现顺序）。
    """
    grouped: dict[tuple, list[Enpt2State]] = {}
    for state in states:
        key = (state.block, *state.sector) if per_block else state.sector
        grouped.setdefault(key, []).append(state)
    return grouped


def _state_rows(states: Iterable[Enpt2State], show_block: bool) -> list[list]:
    rows = []
    for state in states:
        row = [
            None if state.spintwo is None else state.spintwo,
            None if state.symmetry is None else state.symmetry,
            state.state,
            state.ncfg,
            state.ncsf,
            state.evar,
            state.ept,
            state.etot,
        ]
        if show_block:
            row.insert(0, state.block)
        rows.append(row)
    return rows


def print_enpt2_states(
    states: Iterable[Enpt2State],
    tablefmt: str = "grid",
    show_block: bool = True,
) -> None:
    """把态列表按“一行一个态”打印成表格。"""
    states = list(states)
    headers = ["spintwo", "symmetry", "state", "ncfg", "ncsf", "evar", "ept", "etot"]
    if show_block:
        headers.insert(0, "block")
    print(
        tabulate(
            _state_rows(states, show_block),
            headers,
            tablefmt=tablefmt,
            stralign="right",
            numalign="right",
            floatfmt=".8f",
            missingval="-",
        )
    )


# ---------------------------------------------------------------------------
# 命令行入口
# ---------------------------------------------------------------------------


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="提取 iCIPT2 输出中的 iCI_ENPT(2)_NonRela::Info 代码块，按态输出。",
    )
    parser.add_argument("files", nargs="+", help="iCIPT2 输出文件")
    parser.add_argument(
        "--sector",
        action="store_true",
        help="按 (spintwo, symmetry) 扇区分组输出",
    )
    parser.add_argument(
        "--tablefmt",
        default="grid",
        help="tabulate 表格格式（默认 grid，可选 plain/tsv/...）",
    )
    parser.add_argument(
        "--per-block",
        action="store_true",
        help="配合 --sector 时，把不同迭代步的同一扇区分开",
    )
    args = parser.parse_args(argv)

    for path in args.files:
        blocks = read_enpt2_info(path)
        if len(args.files) > 1:
            print(f"=== {path} ===")
        if not blocks:
            print(f"未在 {path} 中找到 {BLOCK_TITLE} 代码块。", file=sys.stderr)
            continue

        states = [state for block in blocks for state in block.states]
        if not args.sector:
            print_enpt2_states(states, tablefmt=args.tablefmt)
            continue

        grouped = group_states_by_sector(states, per_block=args.per_block)
        for key, sub in grouped.items():
            if args.per_block:
                block_no, spintwo, symmetry = key
                print(f"block {block_no}, spintwo {spintwo}, symmetry {symmetry}:")
            else:
                spintwo, symmetry = key
                print(f"spintwo {spintwo}, symmetry {symmetry}:")
            print_enpt2_states(sub, tablefmt=args.tablefmt, show_block=False)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
