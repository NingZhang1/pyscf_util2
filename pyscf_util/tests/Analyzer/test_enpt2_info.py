import os
import tempfile
import unittest

from pyscf_util.Analyzer.iCIPT2.enpt2_info import (
    Enpt2Block,
    Enpt2State,
    extract_enpt2_states,
    group_states_by_sector,
    parse_enpt2_info,
    read_enpt2_info,
)

# 用户给出的多态样例
MULTI_STATE = """\

*************** iCI_ENPT(2)_NonRela::Info ***************

____________________________________________________________________________________________________________________________________
   spintwo |   symmetry |       ncfg |       ncsf |            evar |             ept |            etot 
____________________________________________________________________________________________________________________________________
         0 |          0 |      95932 |     217314 |    -75.45475170 |     -0.00020504 |    -75.45495673 
           |            |            |            |    -75.26317999 |     -0.00024095 |    -75.26342094 
           |            |            |            |    -75.20765095 |     -0.00019754 |    -75.20784849 
         0 |          1 |      52390 |     114358 |    -75.20760436 |     -0.00024546 |    -75.20784982 
____________________________________________________________________________________________________________________________________


*******************************************************
"""

# 旧格式（单态，无 spintwo/symmetry 列），摘自 SO2_test.out
OLD_FORMAT = """\

*************** iCI_ENPT(2)_NonRela::Info ***************

_______________________________________________________________________________________________
      ncfg |       ncsf |            evar |             ept |            etot 
_______________________________________________________________________________________________
    129890 |     231392 |   -547.69018534 |     -0.02555684 |   -547.71574218 
_______________________________________________________________________________________________


*******************************************************
"""


class TestParseMultiState(unittest.TestCase):
    def setUp(self):
        self.blocks = parse_enpt2_info(MULTI_STATE)

    def test_block_count(self):
        self.assertEqual(len(self.blocks), 1)
        self.assertIsInstance(self.blocks[0], Enpt2Block)
        self.assertTrue(self.blocks[0].has_spin_symmetry)

    def test_state_count(self):
        self.assertEqual(len(self.blocks[0].states), 4)
        self.assertTrue(all(isinstance(s, Enpt2State) for s in self.blocks[0]))

    def test_sector_and_state_index(self):
        # 同一 (spintwo, symmetry) 扇区内态编号从 0 重新开始
        got = [(s.spintwo, s.symmetry, s.state) for s in self.blocks[0].states]
        self.assertEqual(got, [(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 1, 0)])

    def test_inherited_columns(self):
        # 第二、三行的 ncfg/ncsf 留空，应继承第一行
        for s in self.blocks[0].states[:3]:
            self.assertEqual(s.ncfg, 95932)
            self.assertEqual(s.ncsf, 217314)

    def test_energies(self):
        first = self.blocks[0].states[0]
        self.assertAlmostEqual(first.evar, -75.45475170, places=10)
        self.assertAlmostEqual(first.ept, -0.00020504, places=10)
        self.assertAlmostEqual(first.etot, -75.45495673, places=10)

        last = self.blocks[0].states[-1]
        self.assertAlmostEqual(last.evar, -75.20760436, places=10)
        self.assertEqual(last.ncfg, 52390)
        self.assertEqual(last.ncsf, 114358)

    def test_line_numbers(self):
        # 行号是 1-based，且随文件递增
        numbers = [s.line_number for s in self.blocks[0].states]
        self.assertEqual(numbers, sorted(numbers))


class TestParseOldFormat(unittest.TestCase):
    def test_single_state(self):
        blocks = parse_enpt2_info(OLD_FORMAT)
        self.assertEqual(len(blocks), 1)
        self.assertFalse(blocks[0].has_spin_symmetry)
        self.assertEqual(len(blocks[0].states), 1)

        state = blocks[0].states[0]
        self.assertIsNone(state.spintwo)
        self.assertIsNone(state.symmetry)
        self.assertEqual(state.state, 0)
        self.assertEqual(state.ncfg, 129890)
        self.assertAlmostEqual(state.etot, -547.71574218, places=10)


class TestMultipleBlocks(unittest.TestCase):
    def test_blocks_are_kept_separate(self):
        text = MULTI_STATE + "\n\n" + OLD_FORMAT
        blocks = parse_enpt2_info(text)
        self.assertEqual(len(blocks), 2)
        self.assertEqual(blocks[0].index, 0)
        self.assertEqual(blocks[1].index, 1)
        self.assertEqual(len(blocks[0].states), 4)
        self.assertEqual(len(blocks[1].states), 1)

        # 每个块的态编号独立
        self.assertEqual([s.state for s in blocks[0].states], [0, 1, 2, 0])
        self.assertEqual([s.state for s in blocks[1].states], [0])

    def test_flat_extract(self):
        text = MULTI_STATE + "\n\n" + MULTI_STATE
        blocks = parse_enpt2_info(text)
        states = [s for b in blocks for s in b.states]
        self.assertEqual(len(states), 8)
        self.assertEqual([s.block for s in states], [0] * 4 + [1] * 4)


class TestRobustness(unittest.TestCase):
    def test_fortran_d_exponent(self):
        text = MULTI_STATE.replace("-0.00020504", "1.23456789D-02")
        blocks = parse_enpt2_info(text)
        self.assertAlmostEqual(blocks[0].states[0].ept, 1.23456789e-2, places=12)

    def test_short_rows_without_spin_columns(self):
        # 某些输出会把 spintwo/symmetry 两列整体省略
        text = MULTI_STATE.replace(
            "           |            |            |            |",
            "              52390 |     114358 |",
        )
        blocks = parse_enpt2_info(text)
        self.assertEqual(len(blocks[0].states), 4)
        # 被省略列时从上一行继承自旋/对称性，ncfg/ncsf 取本行值
        self.assertEqual(
            (blocks[0].states[1].spintwo, blocks[0].states[1].symmetry), (0, 0)
        )
        self.assertEqual(blocks[0].states[1].ncfg, 52390)

    def test_no_block(self):
        self.assertEqual(parse_enpt2_info("nothing here"), [])
        self.assertEqual(parse_enpt2_info(""), [])


class TestGrouping(unittest.TestCase):
    def test_group_by_sector(self):
        blocks = parse_enpt2_info(MULTI_STATE)
        grouped = group_states_by_sector(blocks[0].states)
        self.assertEqual(set(grouped), {(0, 0), (0, 1)})
        self.assertEqual(len(grouped[(0, 0)]), 3)
        self.assertEqual(len(grouped[(0, 1)]), 1)
        self.assertEqual([s.state for s in grouped[(0, 0)]], [0, 1, 2])

    def test_group_per_block(self):
        blocks = parse_enpt2_info(MULTI_STATE + "\n" + MULTI_STATE)
        states = [s for b in blocks for s in b.states]
        grouped = group_states_by_sector(states, per_block=True)
        self.assertEqual(len(grouped), 4)
        self.assertIn((0, 0, 0), grouped)
        self.assertIn((1, 0, 1), grouped)


class TestFileIO(unittest.TestCase):
    def test_read_and_extract(self):
        handle, path = tempfile.mkstemp(suffix=".out", text=True)
        try:
            with os.fdopen(handle, "w") as stream:
                stream.write(MULTI_STATE)
            blocks = read_enpt2_info(path)
            self.assertEqual(len(blocks), 1)
            states = extract_enpt2_states(path)
            self.assertEqual(len(states), 4)
            self.assertEqual(states[0].etot, blocks[0].states[0].etot)
        finally:
            os.remove(path)


if __name__ == "__main__":
    unittest.main()
