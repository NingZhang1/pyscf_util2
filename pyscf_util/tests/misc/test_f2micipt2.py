"""Unit tests for pyscf_util.MeanField.f2micipt2 (no F2M_iCIPT2 binary needed)."""

import os
import tempfile
import unittest

import numpy as np

from pyscf import ao2mo, tools

from pyscf_util.MeanField import f2micipt2


def _determinant_rdm(norb, nocc):
    """Closed-shell determinant: gamma_pq = 2 delta_pq,
    Gamma_pqrs = 4 delta_pq delta_rs - 2 delta_ps delta_qr (PySCF layout)."""
    rdm1 = np.zeros((norb, norb))
    for i in range(nocc):
        rdm1[i, i] = 2.0
    rdm2 = np.zeros((norb,) * 4)
    for p in range(nocc):
        for q in range(nocc):
            for r in range(nocc):
                for s in range(nocc):
                    val = 0.0
                    if p == q and r == s:
                        val += 4.0
                    if p == s and q == r:
                        val -= 2.0
                    rdm2[p, q, r, s] = val
    return rdm1, rdm2


def _write_fake_csv(filename, norb, rdm1, rdm2, with_state=True):
    """Write the RDM csv in the layout produced by the iCIPT2 family codes.

    The value stored under the label (p,q,r,s) is Gamma_{p,q,s,r} in the
    PySCF layout, i.e. rdm2.transpose(0,3,2,1)[p,q,r,s].
    """
    V = rdm2.transpose(0, 3, 2, 1)
    with open(filename, "w") as f:
        if with_state:
            f.write("state,i,j,rdm1\n")
        else:
            f.write("i,j,rdm1\n")
        for i in range(norb):
            for j in range(i + 1):
                if abs(rdm1[i, j]) > 1e-12:
                    if with_state:
                        f.write("0,%d,%d,%#.15g\n" % (i, j, rdm1[i, j]))
                    else:
                        f.write("%d,%d,%#.15g\n" % (i, j, rdm1[i, j]))

    with open(filename.replace("dm1", "dm2"), "w") as f:
        if with_state:
            f.write("state,p,q,r,s,rdm2\n")
        else:
            f.write("p,q,r,s,rdm2\n")
        for p in range(norb):
            for q in range(norb):
                for r in range(norb):
                    for s in range(norb):
                        if abs(V[p, q, r, s]) > 1e-12:
                            if with_state:
                                f.write(
                                    "0,%d,%d,%d,%d,%#.15g\n"
                                    % (p, q, r, s, V[p, q, r, s])
                                )
                            else:
                                f.write(
                                    "%d,%d,%d,%d,%#.15g\n" % (p, q, r, s, V[p, q, r, s])
                                )


class TestRDMReaders(unittest.TestCase):
    def setUp(self):
        self.norb = 4
        self.nocc = 2
        self.nelec = 2 * self.nocc
        self.rdm1, self.rdm2 = _determinant_rdm(self.norb, self.nocc)
        np.random.seed(1)
        h1 = np.random.rand(self.norb, self.norb)
        h1 = h1 + h1.T
        eri = np.random.rand(self.norb, self.norb, self.norb, self.norb)
        # 强制 8 重对称
        eri = ao2mo.restore(1, ao2mo.restore(8, eri, self.norb), self.norb)
        self.h1, self.eri = h1, eri

    def _check(self, with_state):
        with tempfile.TemporaryDirectory() as tmp:
            fname = os.path.join(tmp, "fe2s2_dm1.csv")
            _write_fake_csv(
                fname, self.norb, self.rdm1, self.rdm2, with_state=with_state
            )
            rdm1 = f2micipt2.read_rdm1_csv(fname, self.norb)
            rdm2 = f2micipt2.read_rdm2_csv(fname.replace("dm1", "dm2"), self.norb)

        self.assertTrue(np.allclose(rdm1, self.rdm1, atol=1e-12))
        self.assertTrue(np.allclose(rdm2, self.rdm2, atol=1e-12))

        # 1-RDM 可以由 2-RDM 收缩得到
        rec = np.einsum("pqrr->pq", rdm2) / (self.nelec - 1)
        self.assertTrue(np.allclose(rec, rdm1, atol=1e-10))

        # RDM 能量 = 行列式能量
        e_rdm = f2micipt2.rdm_energy(self.h1, self.eri, rdm1, rdm2, 0.0)
        e_ref = 2.0 * np.einsum("ii->", self.h1[: self.nocc, : self.nocc])
        e_ref += 2.0 * np.einsum(
            "iijj->", self.eri[: self.nocc, : self.nocc, : self.nocc, : self.nocc]
        ) - np.einsum(
            "ijji->", self.eri[: self.nocc, : self.nocc, : self.nocc, : self.nocc]
        )
        self.assertAlmostEqual(e_rdm, e_ref, places=10)

    def test_state_column_format(self):
        self._check(with_state=True)

    def test_state_averaged_format(self):
        self._check(with_state=False)


class TestFCIDUMP(unittest.TestCase):
    def test_write_fcidump_roundtrip(self):
        norb, nelec = 4, 4
        np.random.seed(2)
        h1 = np.random.rand(norb, norb)
        h1 = h1 + h1.T
        eri = np.random.rand(norb, norb, norb, norb)
        # 强制 8 重对称
        eri = ao2mo.restore(1, ao2mo.restore(8, eri, norb), norb)
        ecore = -1.2345

        with tempfile.TemporaryDirectory() as tmp:
            fname = os.path.join(tmp, "FCIDUMP")
            f2micipt2._write_fcidump(fname, h1, eri, norb, nelec, ecore=ecore)
            data = tools.fcidump.read(fname)

        self.assertEqual(data["NORB"], norb)
        self.assertEqual(data["NELEC"], nelec)
        self.assertAlmostEqual(data["ECORE"], ecore, places=12)
        self.assertTrue(np.allclose(data["H1"], h1, atol=1e-12))
        eri1 = ao2mo.restore(1, data["H2"], norb)
        self.assertTrue(np.allclose(eri1, eri, atol=1e-12))


class TestParseEnergies(unittest.TestCase):
    def test_parse_eig(self):
        text = """
*************** Begin MatrixDiag ***************
Init[   0] :   -115.95372006
Eigenvalue :
Eig[   0]  :   -116.47094408
*************** End MatrixDiag ***************
Eig[   0]  :   -116.38405772
"""
        self.assertEqual(f2micipt2.parse_ci_energies(text, 1), [-116.38405772])

    def test_parse_task_final_energies(self):
        text = """
task 0 ene of sparse (init guess) : -115.953720058215
task 0 final energies:
  state  0 : -116.384057717351
  state  1 : -116.200000000000
task 1 final energies:
  state  0 : -116.400000000000
"""
        self.assertEqual(f2micipt2.parse_task_energies(text, 1), [-116.400000000000])
        # 只取最后一个 task 的态
        self.assertEqual(f2micipt2.parse_task_energies(text, 2), [-116.400000000000])
        self.assertIsNone(f2micipt2.parse_task_energies("no energies here"))


class TestOutputParsers(unittest.TestCase):
    LOG = """\
-------------------------------------------
PrimitiveLocalStateLoopUpTable
block_id         : 0
start_orb_indx   : 0
norb             : 1
mx_spintwo       : 20
mx_nopen         : 5
local state info : 
    nelec : 2, spintwo : 0, nopen : 0
        n_local_state                 : 1
        primitive_local_states_offset : 0
    nelec : 1, spintwo : 1, nopen : 1
        n_local_state                 : 1
        primitive_local_states_offset : 1

-------------------------------------------
PrimitiveLocalStateLoopUpTable
block_id         : 1
start_orb_indx   : 1
norb             : 3
local state info : 
    nelec : 6, spintwo : 0, nopen : 0
        n_local_state                 : 1
    nelec : 5, spintwo : 1, nopen : 1
        n_local_state                 : 3

*************** Begin PRINT BLOCK CSF SPACE ***************

spintwo    = 0
nstates    = 1
ncsf       = 10
*********** state  0 ***********
Important block csf:
(   1,  0)(   1,  0)     -0.9 (  90.0%)

*************** End PRINT BLOCK CSF SPACE ***************
"""

    def test_parse_local_states(self):
        blocks = f2micipt2.parse_local_states(self.LOG)
        self.assertEqual(len(blocks), 2)
        self.assertEqual(blocks[0]["block"], 0)
        self.assertEqual(blocks[0]["norb"], 1)
        self.assertEqual(blocks[0]["states"], [(2, 0, 0, 1), (1, 1, 1, 1)])
        self.assertEqual(blocks[1]["block"], 1)
        self.assertEqual(blocks[1]["norb"], 3)
        self.assertEqual(blocks[1]["states"], [(6, 0, 0, 1), (5, 1, 1, 3)])

    def test_parse_important_csfs(self):
        csfs = f2micipt2.parse_important_csfs(self.LOG)
        self.assertEqual(len(csfs), 1)
        self.assertAlmostEqual(csfs[0]["coef"], -0.9)
        self.assertAlmostEqual(csfs[0]["weight"], 90.0)

    def test_extract_important_csf(self):
        text = f2micipt2.extract_important_csf(self.LOG)
        self.assertIn("ncsf       = 10", text)
        self.assertIn("Important block csf", text)
        self.assertNotIn("End PRINT", text)


class TestSolverCLI(unittest.TestCase):
    def test_fe2s2_defaults(self):
        solver = f2micipt2.F2MiCIPT2Solver(executable="/bin/true", workdir="/tmp/x")
        self.assertEqual(solver.norb_model, 20)
        self.assertEqual(solver.mx_spintwo, 10)
        self.assertEqual(solver.mx_ligand_spintwo, 3)
        self.assertEqual(solver.selection, 1)
        self.assertEqual(solver.dumprdm, 2)
        self.assertEqual(solver.dm1_name, "fe2s2_dm1.csv")
        cmd = solver._build_cmd("inp", "FCIDUMP")
        # <exe> <input> <FCIDUMP> <mx_spintwo> <mode> <lig_ext> <lig_s2>
        # <calc_dm> <only_fa> <brute_thr> <hmat_update> <prune>
        self.assertEqual(
            cmd,
            [
                "/bin/true",
                "inp",
                "FCIDUMP",
                "10",
                "1",
                "1",
                "3",
                "1",
                "0",
                "0",
                "1",
                "1",
            ],
        )

    def test_new_options(self):
        solver = f2micipt2.F2MiCIPT2Solver(
            executable="/bin/true",
            workdir="/tmp/x",
            mx_spintwo=6,
            mode=3,
            mx_ligand_ext_level=6,
            mx_ligand_spintwo=2,
            only_ferro_antiferro=True,
            hmat_brute_force_threshold=5000,
            hmat_update=False,
            prune_init_guess=False,
        )
        self.assertEqual(
            solver._build_cmd("inp", "FCIDUMP"),
            [
                "/bin/true",
                "inp",
                "FCIDUMP",
                "6",
                "3",
                "6",
                "2",
                "1",
                "1",
                "5000",
                "0",
                "0",
            ],
        )

    def test_no_rdm(self):
        solver = f2micipt2.F2MiCIPT2Solver(
            executable="/bin/true", workdir="/tmp/x", calc_rdm=False
        )
        self.assertEqual(solver.dumprdm, 0)
        self.assertEqual(solver._build_cmd("inp", "FCIDUMP")[7], "0")
        with self.assertRaises(RuntimeError):
            solver.make_rdm1()

    def test_invalid_arguments(self):
        with self.assertRaises(ValueError):
            f2micipt2.F2MiCIPT2Solver(executable="/bin/true", mode=4)
        with self.assertRaises(ValueError):
            f2micipt2.F2MiCIPT2Solver(executable="/bin/true", selection=2)
        with self.assertRaises(ValueError):
            f2micipt2.F2MiCIPT2Solver(executable="/bin/true", dumprdm=3)
        with self.assertRaises(ValueError):
            # 请求 RDM 但 dumprdm = 0
            f2micipt2.F2MiCIPT2Solver(executable="/bin/true", calc_rdm=True, dumprdm=0)

    def test_input_file(self):
        solver = f2micipt2.F2MiCIPT2Solver(
            executable="/bin/true",
            workdir="/tmp/x",
            cmin="1e-3",
            etol="1e-6",
            mode=2,
            selection=1,
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "f2m.inp")
            solver._write_input(path, 0)
            with open(path) as f:
                text = f.read()
        self.assertIn("selection=1", text)
        self.assertIn("dumprdm=2", text)
        self.assertIn("task=0 0 1 1", text)
        self.assertIn("cmin=1e-3", text)


class TestFe4S4Solver(unittest.TestCase):
    def test_model_defaults(self):
        from pyscf_util.MeanField import f2micipt2_fe4s4 as fe4s4

        solver = fe4s4.F2MiCIPT2Fe4S4Solver(executable="/bin/true", workdir="/tmp/x")
        self.assertEqual(solver.norb_model, 36)
        self.assertEqual(solver.mx_spintwo, 20)
        self.assertEqual(solver.nsegment, "0 0 18 18 0 0")
        self.assertEqual(solver.dm1_name, fe4s4.F2M_DM1_NAME)
        self.assertEqual(solver.dm2_name, fe4s4.F2M_DM2_NAME)

    def test_cli(self):
        from pyscf_util.MeanField import f2micipt2_fe4s4 as fe4s4

        solver = fe4s4.F2MiCIPT2Fe4S4Solver(
            executable="/bin/true",
            workdir="/tmp/x",
            mode=3,
            mx_ligand_ext_level=2,
            mx_ligand_spintwo=1,
            only_ferro_antiferro=True,
        )
        self.assertEqual(
            solver._build_cmd("inp", "FCIDUMP"),
            [
                "/bin/true",
                "inp",
                "FCIDUMP",
                "20",
                "3",
                "2",
                "1",
                "1",
                "1",
                "0",
                "1",
                "1",
            ],
        )

    def test_invalid_mode(self):
        from pyscf_util.MeanField import f2micipt2_fe4s4 as fe4s4

        with self.assertRaises(ValueError):
            fe4s4.F2MiCIPT2Fe4S4Solver(executable="/bin/true", mode=4)


class TestMoGauge(unittest.TestCase):
    @staticmethod
    def _random_orthogonal(n, seed=0):
        rng = np.random.default_rng(seed)
        q, _ = np.linalg.qr(rng.normal(size=(n, n)))
        return q

    def test_max_diag_sq_unitary(self):
        rng = np.random.default_rng(3)
        M = rng.normal(size=(6, 6))
        U = f2micipt2._max_diag_sq_unitary(M)
        self.assertTrue(np.allclose(U.T @ U, np.eye(6), atol=1e-12))

        def objective(x):
            return float(np.sum(np.diag(x) ** 2))

        self.assertGreaterEqual(objective(M @ U) + 1e-12, objective(M))

    def test_gauge_rotate_mo(self):
        n = 10
        mo = self._random_orthogonal(n, seed=5)
        groups = [list(range(5)), list(range(5, 10))]
        mo2 = f2micipt2.gauge_rotate_mo(mo, groups)

        # 仍然是正交矩阵
        self.assertTrue(np.allclose(mo2.T @ mo2, np.eye(n), atol=1e-12))
        # 每个 group 的列空间不变
        for g in groups:
            p1 = mo[:, g] @ mo[:, g].T
            p2 = mo2[:, g] @ mo2[:, g].T
            self.assertTrue(np.allclose(p1, p2, atol=1e-12))
        # 对角权重不下降
        w1, _ = f2micipt2.mo_diagonal_weight(mo)
        w2, _ = f2micipt2.mo_diagonal_weight(mo2)
        self.assertGreaterEqual(w2 + 1e-12, w1)

    def test_gauge_groups(self):
        solver = f2micipt2.F2MiCIPT2Solver(executable="/bin/true", mode=1)
        occ = [i for b in solver.occupied_block_ids for i in solver.model_blocks[b][1]]

        # 默认只旋转双占 (S-S) 组 (严格能量不变)
        groups = solver.gauge_groups()
        self.assertEqual(len(groups), 1)
        self.assertEqual(sorted(groups[0]), sorted(occ))

        # occ_act 额外旋转活性 Fe 组
        occ_act = solver.gauge_groups("occ_act")
        self.assertEqual(len(occ_act), 2)
        self.assertEqual(sorted(occ_act[0]), sorted(occ))
        flat = sorted(i for g in occ_act for i in g)
        self.assertEqual(flat, list(range(20)))
        self.assertEqual(sorted(solver.gauge_groups("occ")[0]), sorted(occ))

        # mode 3 (LMCT): 全部轨道都是活性轨道, 不做 occ/act gauge
        solver.mode = 3
        self.assertIsNone(solver.gauge_groups())

        # 每个局域块一组总是存在
        local = solver.gauge_groups("local")
        self.assertEqual(len(local), len(solver.model_blocks))
        self.assertEqual(sorted(i for g in local for i in g), list(range(20)))

    def test_mode3_gauge_mo_unchanged(self):
        solver = f2micipt2.F2MiCIPT2Solver(executable="/bin/true", mode=3)
        mo = self._random_orthogonal(20, seed=7)
        mo2 = solver.gauge_mo_coeff(mo)
        self.assertTrue(np.allclose(mo, mo2))

    def test_fe4s4_gauge_groups(self):
        from pyscf_util.MeanField import f2micipt2_fe4s4 as fe4s4

        solver = fe4s4.F2MiCIPT2Fe4S4Solver(executable="/bin/true", mode=1)
        groups = solver.gauge_groups()
        self.assertEqual(len(groups), 1)  # 默认只旋转双占 S 组
        groups2 = solver.gauge_groups("occ_act")
        self.assertEqual(len(groups2), 2)
        flat = sorted(i for g in groups2 for i in g)
        self.assertEqual(flat, list(range(36)))
        self.assertEqual(len(set(flat)), 36)
        self.assertEqual(
            sorted(solver.occupied_block_ids + solver.active_block_ids),
            list(range(len(solver.model_blocks))),
        )


if __name__ == "__main__":
    unittest.main()
