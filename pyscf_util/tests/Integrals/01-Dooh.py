import unittest
import numpy as np
from pyscf import gto, scf, ao2mo
from pyscf.tools import fcidump
from functools import reduce
import os

from pyscf_util.misc.mole import get_orbsym
from pyscf_util.Integrals.integral_Dooh import (
    get_symmetry_adapted_basis_Dooh,
    FCIDUMP_Dooh,
)


class TestC2Molecule(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mol = gto.M(
            verbose=0,
            atom="""
            C   0.000000000000       0.000000000000      -0.621265
            C   0.000000000000       0.000000000000       0.621265
            """,
            basis={"C": "cc-pvdz"},
            spin=0,
            charge=0,
            symmetry="dooh",
        )
        cls.mol.build()
        cls.mf = scf.RHF(cls.mol)
        cls.mf.kernel()

    def test_scf_convergence(self):
        self.assertTrue(self.mf.converged)

    def test_symmetry(self):
        orbsym_ID, orbsym = get_orbsym(self.mol, self.mf.mo_coeff)
        self.assertIsNotNone(orbsym_ID)
        self.assertIsNotNone(orbsym)

    def test_symmetry_adapted_basis(self):
        basis_trans, Lz, parity = get_symmetry_adapted_basis_Dooh(
            self.mol, self.mf.mo_coeff
        )
        self.assertIsNotNone(basis_trans)
        self.assertIsNotNone(Lz)
        self.assertIsNotNone(parity)

    def test_h1e_symmetry(self):
        basis_trans, Lz, parity = get_symmetry_adapted_basis_Dooh(
            self.mol, self.mf.mo_coeff
        )
        h1e = reduce(
            np.dot, (self.mf.mo_coeff.T, self.mf.get_hcore(), self.mf.mo_coeff)
        )
        h1e_adapted = reduce(np.dot, (basis_trans.H, h1e, basis_trans))

        for i in range(self.mol.nao):
            for j in range(self.mol.nao):
                if abs(h1e_adapted[i, j]) > 1e-10:
                    self.assertAlmostEqual(Lz[i], Lz[j], places=10)
                    self.assertEqual(parity[i], parity[j])

    def test_int2e_symmetry(self):
        basis_trans, Lz, parity = get_symmetry_adapted_basis_Dooh(
            self.mol, self.mf.mo_coeff
        )
        int2e_full = ao2mo.full(self.mol, self.mf.mo_coeff, aosym="1").reshape(
            (self.mol.nao,) * 4
        )

        int2e_full = np.einsum("ijkl,ip->pjkl", int2e_full, basis_trans.conj())
        int2e_full = np.einsum("pjkl,jq->pqkl", int2e_full, basis_trans)
        int2e_full = np.einsum("pqkl,kr->pqrl", int2e_full, basis_trans.conj())
        int2e_full = np.einsum("pqrl,ls->pqrs", int2e_full, basis_trans)

        for i in range(self.mol.nao):
            for j in range(self.mol.nao):
                for k in range(self.mol.nao):
                    for l in range(self.mol.nao):
                        if abs(int2e_full[i, j, k, l]) > 1e-10:
                            self.assertLess(abs(int2e_full[i, j, k, l].imag), 1e-10)
                            self.assertLess(abs(Lz[i] + Lz[k] - Lz[j] - Lz[l]), 1e-10)
                            self.assertEqual(
                                parity[i] * parity[j] * parity[k] * parity[l], 1
                            )

    def test_fcidump_creation(self):
        FCIDUMP_Dooh(self.mol, self.mf, "FCIDUMP_C2_test")
        # Here you might want to add checks to ensure the FCIDUMP file was created correctly

    def test_fcidump_full_creation(self):
        filename = "FCIDUMP_C2_FULL_test"
        nmo = self.mf.mo_coeff.shape[1]
        nelec = self.mol.nelectron
        ms = 0
        tol = 1e-10
        nuc = self.mol.get_enuc()

        orbsym_ID, _ = get_orbsym(self.mol, self.mf.mo_coeff)
        basis_trans, _, _ = get_symmetry_adapted_basis_Dooh(self.mol, self.mf.mo_coeff)

        h1e = reduce(
            np.dot, (self.mf.mo_coeff.T, self.mf.get_hcore(), self.mf.mo_coeff)
        )
        h1e_adapted = reduce(np.dot, (basis_trans.H, h1e, basis_trans)).real

        int2e_full = ao2mo.full(self.mol, self.mf.mo_coeff, aosym="1").reshape(
            (self.mol.nao,) * 4
        )
        int2e_full = np.einsum("ijkl,ip->pjkl", int2e_full, basis_trans.conj())
        int2e_full = np.einsum("pjkl,jq->pqkl", int2e_full, basis_trans)
        int2e_full = np.einsum("pqkl,kr->pqrl", int2e_full, basis_trans.conj())
        int2e_full = np.einsum("pqrl,ls->pqrs", int2e_full, basis_trans).real

        with open(filename, "w") as fout:
            fcidump.write_head(fout, nmo, nelec, ms, orbsym_ID)
            fcidump.write_eri(fout, int2e_full, nmo, tol=tol)
            fcidump.write_hcore(fout, h1e_adapted, nmo, tol=tol)
            fout.write(f"{nuc:.16f}  0  0  0  0\n")

        self.assertTrue(os.path.exists(filename))


if __name__ == "__main__":
    unittest.main()
