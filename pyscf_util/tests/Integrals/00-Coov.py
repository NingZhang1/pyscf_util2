import unittest
import numpy as np
from pyscf import gto, scf, ao2mo
from functools import reduce

from pyscf_util.Integrals.integral_Coov import (
    get_symmetry_adapted_basis_Coov,
    FCIDUMP_Coov,
)
from pyscf_util.misc.mole import get_orbsym


class TestCOMolecule(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mol = gto.M(
            verbose=0,
            atom="""
            C   0.000000000000       0.000000000000      -0.621265
            O   0.000000000000       0.000000000000       0.621265
            """,
            basis={"C": "cc-pvdz", "O": "cc-pvdz"},
            spin=0,
            charge=0,
            symmetry="coov",
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
        basis_trans, Lz = get_symmetry_adapted_basis_Coov(self.mol, self.mf.mo_coeff)
        self.assertIsNotNone(basis_trans)
        self.assertIsNotNone(Lz)

    def test_h1e_symmetry(self):
        basis_trans, Lz = get_symmetry_adapted_basis_Coov(self.mol, self.mf.mo_coeff)
        h1e = reduce(
            np.dot, (self.mf.mo_coeff.T, self.mf.get_hcore(), self.mf.mo_coeff)
        )
        h1e_adapted = reduce(np.dot, (basis_trans.H, h1e, basis_trans))

        for i in range(self.mol.nao):
            for j in range(self.mol.nao):
                if abs(h1e_adapted[i, j]) > 1e-10:
                    self.assertAlmostEqual(Lz[i], Lz[j], places=10)

    def test_int2e_symmetry(self):
        basis_trans, Lz = get_symmetry_adapted_basis_Coov(self.mol, self.mf.mo_coeff)
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

    def test_fcidump_creation(self):
        FCIDUMP_Coov(self.mol, self.mf, "FCIDUMP_CO_test")
        # Here you might want to add checks to ensure the FCIDUMP file was created correctly


if __name__ == "__main__":
    unittest.main()
