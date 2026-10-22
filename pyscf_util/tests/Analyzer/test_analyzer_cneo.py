import os
import tempfile
import unittest

from pyscf_util.Analyzer.iCIPT2.analyzer_cneo import (
    CneoBlock,
    CneoRoot,
    CneoRow,
    CneoSpace,
    CneoState,
    cneo_statistics,
    extract_cneo_states,
    iter_cneo_states,
    parse_cneo_icipt2,
    read_cneo_icipt2,
    split_by_kind,
)

# 用户给出的实际样例（ES + PES，各一个 Space、一个态）
REAL_SAMPLE = """\
*******************************************************
*       cNEO-iCIPT2 ES       *
-------------------------------------------------------
Space : spintwo 0 irrep 0
--------------------------------------------------------------------------------------------------------------
|       cmin |       ncsf |       ncfg | (                eiCI,                ePT2,                etot) |
--------------------------------------------------------------------------------------------------------------
| 1.0000e-04 |      10908 |       7823 | (    -99.247548998579,     -0.010738803995,    -99.258287802574) |
| 5.0000e-05 |      27444 |      18772 | (    -99.250658681049,     -0.007668211512,    -99.258326892561) |
| 3.0000e-05 |      58403 |      38022 | (    -99.252792949383,     -0.005566410532,    -99.258359359915) |
| 1.5000e-05 |     160098 |      97397 | (    -99.255187847603,     -0.003205463762,    -99.258393311365) |
| 9.0000e-06 |     310755 |     178306 | (    -99.256382869736,     -0.002033695199,    -99.258416564935) |
| 7.0000e-06 |     419589 |     234857 | (    -99.256810535424,     -0.001615301969,    -99.258425837393) |
-------------------------------------------------------
Nuclear Kinetic Energy
State 0
1.0000e-04 |       0.468578994258 
5.0000e-05 |       0.467499836025 
3.0000e-05 |       0.467001347937 
1.5000e-05 |       0.466585790658 
9.0000e-06 |       0.466421701795 
7.0000e-06 |       0.466365059907 
-------------------------------------------------------
*       cNEO-iCIPT2 PES      *
-------------------------------------------------------
Space : spintwo 0 irrep 0
--------------------------------------------------------------------------------------------------------------
|       cmin |       ncsf |       ncfg | (                eiCI,                ePT2,                etot) |
--------------------------------------------------------------------------------------------------------------
| 1.0000e-04 |      10908 |       7823 | (    -99.716127992837,     -0.010738803995,    -99.726866796832) |
| 5.0000e-05 |      27444 |      18772 | (    -99.718158517074,     -0.007668211512,    -99.725826728586) |
| 3.0000e-05 |      58403 |      38022 | (    -99.719794297319,     -0.005566410532,    -99.725360707851) |
| 1.5000e-05 |     160098 |      97397 | (    -99.721773638261,     -0.003205463762,    -99.724979102024) |
| 9.0000e-06 |     310755 |     178306 | (    -99.722804571531,     -0.002033695199,    -99.724838266731) |
| 7.0000e-06 |     419589 |     234857 | (    -99.723175595331,     -0.001615301969,    -99.724790897300) |
-------------------------------------------------------
--------------------------------------------------------------------------------------------------------------
*******************************************************
"""

# 多 Space、多根、多核动能态的样例
MULTI = """\
*******************************************************
*       cNEO-iCIPT2 ES       *
-------------------------------------------------------
Space : spintwo 0 irrep 0
-------------------------------------------------------
|       cmin |       ncsf |       ncfg | (   eiCI,   ePT2,   etot) |
-------------------------------------------------------
| 1.0000e-04 |         10 |          5 | ( 1.0, 0.1, 1.1 ) ( 2.0, 0.2, 2.2 ) |
| 5.0000e-05 |         20 |          9 | ( 1.5, 0.05, 1.55 ) ( 2.5, 0.25, 2.75 ) |
-------------------------------------------------------
Nuclear Kinetic Energy
State 0
1.0000e-04 |       0.50
5.0000e-05 |       0.40
State 1
1.0000e-04 |       0.60
5.0000e-05 |       0.70
-------------------------------------------------------
Space : spintwo 0 irrep 1
-------------------------------------------------------
|       cmin |       ncsf |       ncfg | (   eiCI,   ePT2,   etot) |
-------------------------------------------------------
| 1.0000e-04 |         11 |          6 | ( 3.0, 0.3, 3.3 ) |
-------------------------------------------------------
Nuclear Kinetic Energy
State 0
1.0000e-04 |       0.90
-------------------------------------------------------
*       cNEO-iCIPT2 PES      *
-------------------------------------------------------
Space : spintwo 0 irrep 0
-------------------------------------------------------
|       cmin |       ncsf |       ncfg | (   eiCI,   ePT2,   etot) |
-------------------------------------------------------
| 1.0000e-04 |         10 |          5 | ( 4.0, 0.4, 4.4 ) |
-------------------------------------------------------
-------------------------------------------------------
*******************************************************
"""


class TestRealSample(unittest.TestCase):
    def setUp(self):
        self.blocks = parse_cneo_icipt2(REAL_SAMPLE)

    def test_blocks(self):
        self.assertEqual(len(self.blocks), 2)
        self.assertTrue(all(isinstance(b, CneoBlock) for b in self.blocks))
        self.assertEqual([b.kind for b in self.blocks], ["ES", "PES"])
        self.assertEqual([b.index for b in self.blocks], [0, 1])

    def test_spaces(self):
        self.assertEqual(len(self.blocks[0].spaces), 1)
        space = self.blocks[0].spaces[0]
        self.assertIsInstance(space, CneoSpace)
        self.assertEqual((space.spintwo, space.irrep), (0, 0))
        self.assertEqual(space.kind, "ES")
        self.assertEqual(space.block, 0)

    def test_rows(self):
        space = self.blocks[0].spaces[0]
        self.assertEqual(len(space.rows), 6)
        self.assertTrue(all(isinstance(r, CneoRow) for r in space.rows))
        self.assertEqual(len(space.rows[0].roots), 1)
        self.assertIsInstance(space.rows[0].roots[0], CneoRoot)
        self.assertEqual(space.rows[0].ncsf, 10908)
        self.assertEqual(space.rows[0].ncfg, 7823)

    def test_first_and_last_values(self):
        first = self.blocks[0].spaces[0].rows[0]
        self.assertAlmostEqual(first.cmin, 1.0e-4, places=16)
        self.assertAlmostEqual(first.roots[0].eiCI, -99.247548998579, places=10)
        self.assertAlmostEqual(first.roots[0].ePT2, -0.010738803995, places=12)
        self.assertAlmostEqual(first.roots[0].etot, -99.258287802574, places=10)

        pes = self.blocks[1].spaces[0].rows[-1]
        self.assertAlmostEqual(pes.roots[0].etot, -99.724790897300, places=10)

    def test_nuclear_kinetic(self):
        space = self.blocks[0].spaces[0]
        self.assertIn(0, space.nuclear_kinetic)
        self.assertEqual(len(space.nuclear_kinetic[0]), 6)
        cmin, value = space.nuclear_kinetic[0][0]
        self.assertAlmostEqual(cmin, 1.0e-4, places=16)
        self.assertAlmostEqual(value, 0.468578994258, places=12)
        self.assertEqual(space.nuclear_kinetic_map(0)[1.0e-4], value)
        # PES 没有核动能段
        self.assertEqual(self.blocks[1].spaces[0].nuclear_kinetic, {})

    def test_statistics(self):
        stats = cneo_statistics(self.blocks)
        self.assertEqual(stats["ES"]["blocks"], 1)
        self.assertEqual(stats["ES"]["spaces"], 1)
        self.assertEqual(stats["ES"]["cmin_points"], 6)
        self.assertEqual(stats["ES"]["states"], 1)
        self.assertEqual(stats["ES"]["records"], 6)
        self.assertEqual(stats["PES"]["records"], 6)

    def test_split_by_kind(self):
        grouped = split_by_kind(self.blocks)
        self.assertEqual(list(grouped), ["ES", "PES"])
        self.assertEqual(len(grouped["ES"]), 1)
        self.assertEqual(len(grouped["PES"]), 1)


class TestMultiSpaceMultiRoot(unittest.TestCase):
    def setUp(self):
        self.blocks = parse_cneo_icipt2(MULTI)
        self.states = list(iter_cneo_states(self.blocks))

    def test_block_and_space_count(self):
        self.assertEqual(len(self.blocks), 2)
        self.assertEqual(len(self.blocks[0].spaces), 2)
        self.assertEqual(len(self.blocks[1].spaces), 1)

    def test_multiple_roots_per_row(self):
        first_space = self.blocks[0].spaces[0]
        self.assertEqual(len(first_space.rows[0].roots), 2)
        self.assertEqual(first_space.nstates, 2)
        self.assertAlmostEqual(first_space.rows[1].roots[1].eiCI, 2.5, places=12)

    def test_flattened_states(self):
        # ES: 2 个 space，第一 space 2 cmin × 2 态 = 4，第二 space 1 cmin × 1 态 = 1
        # PES: 1 space，1 cmin × 1 态 = 1
        self.assertEqual(len(self.states), 6)
        self.assertEqual(sum(1 for s in self.states if s.kind == "ES"), 5)
        self.assertEqual(sum(1 for s in self.states if s.kind == "PES"), 1)

    def test_state_numbering(self):
        first_two = self.states[:2]
        self.assertEqual([s.state for s in first_two], [0, 1])
        self.assertEqual(first_two[0].cmin, first_two[1].cmin)
        self.assertAlmostEqual(first_two[1].etot, 2.2, places=12)

    def test_nuclear_kinetic_mapping(self):
        # state 0 / state 1 分别对应各自 cmin 的核动能
        self.assertAlmostEqual(self.states[0].nuclear_kinetic, 0.50, places=12)
        self.assertAlmostEqual(self.states[1].nuclear_kinetic, 0.60, places=12)
        self.assertAlmostEqual(self.states[2].nuclear_kinetic, 0.40, places=12)
        self.assertAlmostEqual(self.states[3].nuclear_kinetic, 0.70, places=12)
        # 第二个 space（irrep 1）
        self.assertAlmostEqual(self.states[4].nuclear_kinetic, 0.90, places=12)
        # PES 没有核动能
        self.assertIsNone(self.states[5].nuclear_kinetic)

    def test_state_fields(self):
        state = self.states[0]
        self.assertIsInstance(state, CneoState)
        self.assertEqual(
            (state.kind, state.block, state.spintwo, state.irrep),
            ("ES", 0, 0, 0),
        )
        self.assertEqual((state.state, state.ncsf, state.ncfg), (0, 10, 5))

    def test_statistics(self):
        stats = cneo_statistics(self.blocks)
        # ES: 2 个 space；state 数 2 + 1 = 3；cmin 行 2 + 1 = 3；记录 4 + 1 = 5
        self.assertEqual(stats["ES"]["spaces"], 2)
        self.assertEqual(stats["ES"]["states"], 3)
        self.assertEqual(stats["ES"]["cmin_points"], 3)
        self.assertEqual(stats["ES"]["records"], 5)
        self.assertEqual(stats["PES"]["states"], 1)
        self.assertEqual(stats["PES"]["records"], 1)


class TestExtractAndFileIO(unittest.TestCase):
    def _write(self, text):
        handle, path = tempfile.mkstemp(suffix=".out", text=True)
        with os.fdopen(handle, "w") as stream:
            stream.write(text)
        return path

    def test_extract_by_kind(self):
        path = self._write(REAL_SAMPLE)
        try:
            es = extract_cneo_states(path, kind="ES")
            pes = extract_cneo_states(path, kind="PES")
            both = extract_cneo_states(path)
            self.assertEqual(len(es), 6)
            self.assertEqual(len(pes), 6)
            self.assertEqual(len(both), 12)
            self.assertTrue(all(s.kind == "ES" for s in es))
            self.assertTrue(all(s.kind == "PES" for s in pes))
        finally:
            os.remove(path)

    def test_read_file(self):
        path = self._write(REAL_SAMPLE)
        try:
            blocks = read_cneo_icipt2(path)
            self.assertEqual(len(blocks), 2)
            self.assertEqual(blocks[0].kind, "ES")
        finally:
            os.remove(path)

    def test_no_block(self):
        self.assertEqual(parse_cneo_icipt2("nothing here"), [])
        self.assertEqual(parse_cneo_icipt2(""), [])


if __name__ == "__main__":
    unittest.main()
