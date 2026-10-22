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
                                    "%d,%d,%d,%d,%#.15g\n"
                                    % (p, q, r, s, V[p, q, r, s])
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
            rdm2 = f2micipt2.read_rdm2_csv(
                fname.replace("dm1", "dm2"), self.norb
            )

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
    def test_parse(self):
        text = """
*************** Begin MatrixDiag ***************
Init[   0] :   -115.95372006
Eigenvalue :
Eig[   0]  :   -116.47094408
*************** End MatrixDiag ***************
Eig[   0]  :   -116.38405772
"""
        self.assertEqual(f2micipt2.parse_ci_energies(text, 1), [-116.38405772])


if __name__ == "__main__":
    unittest.main()
