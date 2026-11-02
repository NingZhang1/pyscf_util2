#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for the optimized four-component integral export.

Run with::

    python pyscf_util/tests/Relativisitc/00-4C-Integrals.py

Set ``RUN_SLOW=1`` to include the cc-pVDZ-DK performance comparison.

The old implementation (:mod:`pyscf_util.Relativisitc.integral_4C`) is used
as the oracle.  All checks compare the *reconstructed* FCIDUMP content, not
the raw bytes, because the AO->MO transformation path is different (blocked,
AO-symmetry aware transforms instead of plain ``einsum``) and therefore
differs in the last floating-point digits.
"""

import os
import shutil
import tempfile
import time
import unittest

import numpy as np
import pyscf
from pyscf import gto, lib, scf

from pyscf_util.Relativisitc import integral_4C as old_mod
from pyscf_util.Relativisitc import integral_4C_opt as opt_mod

RUN_SLOW = os.environ.get("RUN_SLOW", "0") not in ("0", "", "False", "false")

WATER = "O 0 0 0; H 0 -0.757 0.587; H 0 0.757 0.587"


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def parse_fcidump(path):
    """Return ``(records, h1e, ecore, orbsym, nelec, norb)``.

    ``records`` maps ``(i, j, k, l)`` (1-based, as written) to a complex
    value.  Coulomb and Breit records are distinguished by their index range,
    exactly as in the file.
    """
    import re

    records = {}
    h1e = {}
    ecore = None
    orbsym = None
    nelec = None
    norb = None
    with open(path) as fin:
        lines = fin.read().splitlines()
    end = 0
    for end, line in enumerate(lines):
        if "&END" in line:
            break
    head = "\n".join(lines[: end + 1])
    m = re.search(r"NORB\s*=\s*(\d+)", head)
    if m:
        norb = int(m.group(1))
    m = re.search(r"NELEC\s*=\s*(\d+)", head)
    if m:
        nelec = int(m.group(1))
    m = re.search(r"ORBSYM\s*=\s*([0-9,\s]+)", head)
    if m:
        orbsym = [int(x) for x in re.split(r"[,\s]+", m.group(1).strip()) if x]
    for line in lines[end + 1 :]:
        t = line.split()
        if len(t) < 5:
            continue
        if len(t) == 6:
            value = complex(float(t[0]), float(t[1]))
            i, j, k, l = (int(x) for x in t[2:6])
        else:
            value = complex(float(t[0]), 0.0)
            i, j, k, l = (int(x) for x in t[1:5])
        if k == 0:
            if j == 0:
                ecore = value
            else:
                h1e[(i, j)] = value
            continue
        records[(i, j, k, l)] = value
    return records, h1e, ecore, orbsym, nelec, norb


def _same_records(test, ref, got, tol=1e-10, label=""):
    test.assertEqual(set(ref), set(got), label)
    if not ref:
        return
    delta = max(abs(ref[k] - got[k]) for k in ref)
    test.assertLessEqual(delta, tol, label)


class _FourComponentCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp(prefix="rela4c_opt_test_")
        cls.mol = gto.M(
            atom=WATER,
            basis="3-21g",
            verbose=0,
            charge=0,
            spin=0,
            symmetry="c2v",
        )
        cls.mol.build()
        cls.mf = scf.dhf.RDHF(cls.mol)
        cls.mf.conv_tol = 1e-10
        cls.mf.kernel()
        cls.mf.with_breit = True
        cls.npes = 20

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def path(self, name):
        return os.path.join(self.tmpdir, name)


# ---------------------------------------------------------------------------
# 1. dump equivalence
# ---------------------------------------------------------------------------


class TestDumpEquivalence(_FourComponentCase):
    def _run_pair(self, mode, with_breit, is_complex, **kwargs):
        ref = self.path(
            "ref_%s_%d_%d_%s.txt" % (mode, with_breit, is_complex, hash(str(kwargs)))
        )
        new = self.path(
            "new_%s_%d_%d_%s.txt" % (mode, with_breit, is_complex, hash(str(kwargs)))
        )
        old_mod.FCIDUMP_Rela4C(
            self.mol,
            self.mf,
            with_breit=with_breit,
            npes=self.npes,
            filename=ref,
            mode=mode,
            IsComplex=is_complex,
            **kwargs,
        )
        opt_mod.FCIDUMP_Rela4C_opt(
            self.mol,
            self.mf,
            with_breit=with_breit,
            npes=self.npes,
            filename=new,
            mode=mode,
            IsComplex=is_complex,
            **kwargs,
        )
        return parse_fcidump(ref), parse_fcidump(new)

    def test_incore(self):
        for with_breit in (False, True):
            for is_complex in (False, True):
                ref, new = self._run_pair("incore", with_breit, is_complex)
                self._compare(ref, new, "incore breit=%s complex=%s" % (with_breit, is_complex))

    def test_outcore(self):
        for with_breit in (False, True):
            for is_complex in (False, True):
                ref, new = self._run_pair("outcore", with_breit, is_complex)
                self._compare(ref, new, "outcore breit=%s complex=%s" % (with_breit, is_complex))

    def test_incore_equals_outcore(self):
        for with_breit in (False, True):
            _, inc = self._run_pair("incore", with_breit, True)
            _, out = self._run_pair("outcore", with_breit, True)
            self._compare(inc, out, "incore vs outcore breit=%s" % with_breit)

    def test_debug_full_dump(self):
        ref, new = self._run_pair("incore", True, True, debug=True)
        self._compare(ref, new, "debug s1 dump")
        self.assertEqual(len(ref[0]), len(new[0]))
        # the s1 dump must be the full (un-symmetrized) set
        self.assertGreater(len(new[0]), len(self._run_pair("incore", True, True)[1][0]))

    def test_no_2e(self):
        ref, new = self._run_pair("incore", True, True, no_2e=True)
        self.assertEqual(len(ref[0]), 0)
        self.assertEqual(set(ref[1]), set(new[1]))
        self.assertAlmostEqual(ref[2], new[2], places=12)

    def test_paldus_filters_incore(self):
        for kwargs in (
            {"retain_paldus_123": False},
            {"retain_paldus_456": False},
            {"retain_paldus_123": False, "retain_paldus_456": False},
            {"only_dg": True},
        ):
            ref, new = self._run_pair("incore", True, True, **kwargs)
            self._compare(ref, new, "incore %s" % kwargs)

    def test_paldus_filters_outcore_consistent_with_incore(self):
        # NOTE: the old *outcore* dump silently ignored the paldus/dg filters
        #       (a bug); the optimized outcore honors them, matching the
        #       in-core path.  The default (all filters off) is identical.
        for kwargs in (
            {"retain_paldus_123": False},
            {"only_dg": True},
        ):
            _, inc = self._run_pair("incore", True, True, **kwargs)
            _, out = self._run_pair("outcore", True, True, **kwargs)
            self._compare(inc, out, "outcore filters %s" % kwargs)

    def test_with_nes(self):
        ref = self.path("ref_nes.txt")
        new = self.path("new_nes.txt")
        old_mod.FCIDUMP_Rela4C(
            self.mol, self.mf, with_breit=False, npes=10, filename=ref,
            mode="incore", IsComplex=True, with_nes=True,
        )
        opt_mod.FCIDUMP_Rela4C_opt(
            self.mol, self.mf, with_breit=False, npes=10, filename=new,
            mode="incore", IsComplex=True, with_nes=True,
        )
        self._compare(parse_fcidump(ref), parse_fcidump(new), "with_nes")

    def _compare(self, ref, new, label):
        rec_ref, h1_ref, e_ref, sym_ref, n_ref, norb_ref = ref
        rec_new, h1_new, e_new, sym_new, n_new, norb_new = new
        _same_records(self, rec_ref, rec_new, 1e-10, label)
        self.assertEqual(set(h1_ref), set(h1_new), label)
        if h1_ref:
            self.assertLessEqual(
                max(abs(h1_ref[k] - h1_new[k]) for k in h1_ref), 1e-12, label
            )
        self.assertAlmostEqual(e_ref, e_new, places=12)
        self.assertEqual(sym_ref, sym_new, label)
        self.assertEqual(n_ref, n_new, label)
        self.assertEqual(norb_ref, norb_new, label)


# ---------------------------------------------------------------------------
# 2. tensor-level checks of the symmetries that the optimization relies on
# ---------------------------------------------------------------------------


class TestSymmetryAssumptions(_FourComponentCase):
    def _spinor_blocks(self, npes=12):
        n2c = self.mol.nao_2c()
        mo = np.asarray(self.mf.mo_coeff)[:, n2c : n2c + npes]
        return mo[:n2c, :], mo[n2c:, :]

    def _ao_mo(self, intor, coeffs, aosym="s1"):
        from pyscf.ao2mo import r_outcore

        n = coeffs[0].shape[1]
        path = self.path("sym_%s_%s_%d.h5" % (intor, aosym, id(coeffs)))
        r_outcore.general(self.mol, coeffs, path, intor=intor, aosym=aosym, verbose=0)
        import h5py

        with h5py.File(path, "r") as h5:
            arr = np.array(h5["eri_mo"]).reshape(n, n, n, n)
        os.remove(path)
        return arr

    def test_aosym_s4_coulomb(self):
        """4-fold AO symmetry must be exact for the Coulomb intors."""
        L, S = self._spinor_blocks()
        for intor, coeffs in (
            ("int2e_spinor", (L, L, L, L)),
            ("int2e_spsp1spsp2_spinor", (S, S, S, S)),
            ("int2e_spsp2_spinor", (L, L, S, S)),
            ("int2e_spsp1_spinor", (S, S, L, L)),
        ):
            ref = self._ao_mo(intor, coeffs, "s1")
            got = self._ao_mo(intor, coeffs, "s4")
            self.assertLessEqual(np.abs(ref - got).max(), 1e-11, intor)

    def test_aosym_s4_invalid_for_breit(self):
        """The Breit/Gaunt intors must NOT use s4 (documented finding)."""
        L, S = self._spinor_blocks()
        intor = "int2e_breit_ssp1ssp2_spinor"
        ref = self._ao_mo(intor, (L, S, L, S), "s1")
        got = self._ao_mo(intor, (L, S, L, S), "s4")
        scale = max(1e-30, np.abs(ref).max())
        self.assertGreater(np.abs(ref - got).max() / scale, 1e-3, intor)

    def test_breit_term_relations(self):
        """t2 = conj(t1[j,i,l,k]) and t4 = conj(t3[j,i,l,k])."""
        L, S = self._spinor_blocks(npes=10)
        c1 = 0.5 / lib.param.LIGHT_SPEED

        def transform(ao, a, b, c, d):
            return lib.einsum(
                "pqrs,pi,qj,rk,sl->ijkl", ao, a.conj(), b, c.conj(), d
            )

        for intor_lsls, intor_lssl in (
            ("int2e_breit_ssp1ssp2_spinor", "int2e_breit_ssp1sps2_spinor"),
            ("int2e_ssp1ssp2_spinor", "int2e_ssp1sps2_spinor"),
        ):
            ao1 = self.mol.intor(intor_lsls) * c1**2
            ao3 = self.mol.intor(intor_lssl) * c1**2
            t1 = transform(ao1, L, S, L, S)
            t2 = transform(ao1.conj().transpose(1, 0, 3, 2), S, L, S, L)
            t3 = transform(ao3, L, S, S, L)
            t4 = transform(ao3.transpose(2, 3, 0, 1), S, L, L, S)
            self.assertLessEqual(
                np.abs(t1 - t2.transpose(1, 0, 3, 2).conj()).max(), 1e-12
            )
            self.assertLessEqual(
                np.abs(t3 - t4.transpose(1, 0, 3, 2).conj()).max(), 1e-12
            )

    def test_incore_tensor_matches_old(self):
        ref_c, ref_b = old_mod.FCIDUMP_Rela4C(
            self.mol, self.mf, with_breit=True, npes=self.npes,
            filename=self.path("t_ref.txt"), mode="incore", IsComplex=True,
            debug=True, no_2e=True,
        )
        new_c, new_b = opt_mod.FCIDUMP_Rela4C_opt(
            self.mol, self.mf, with_breit=True, npes=self.npes,
            filename=self.path("t_new.txt"), mode="incore", IsComplex=True,
            debug=True, no_2e=True,
        )
        self.assertLessEqual(np.abs(ref_c - new_c).max(), 1e-11)
        self.assertLessEqual(np.abs(ref_b - new_b).max(), 1e-11)


# ---------------------------------------------------------------------------
# 3. optional performance check
# ---------------------------------------------------------------------------


@unittest.skipUnless(RUN_SLOW, "set RUN_SLOW=1 to run the performance check")
class TestPerformance(unittest.TestCase):
    def test_large_case(self):
        tmp = tempfile.mkdtemp(prefix="rela4c_perf_")
        try:
            mol = gto.M(
                atom=WATER, basis="cc-pvdz-dk", verbose=0,
                charge=0, spin=0, symmetry="c2v",
            )
            mol.build()
            mf = scf.dhf.RDHF(mol)
            mf.conv_tol = 1e-9
            mf.kernel()
            mf.with_breit = True
            npes = min(30, mol.nao_2c())
            report = []
            for with_breit in (False, True):
                for mode in ("incore", "outcore"):
                    ref = os.path.join(tmp, "o_%s_%d.txt" % (mode, with_breit))
                    new = os.path.join(tmp, "n_%s_%d.txt" % (mode, with_breit))
                    t = time.time()
                    old_mod.FCIDUMP_Rela4C(
                        mol, mf, with_breit=with_breit, npes=npes,
                        filename=ref, mode=mode, IsComplex=True,
                    )
                    t_old = time.time() - t
                    t = time.time()
                    opt_mod.FCIDUMP_Rela4C_opt(
                        mol, mf, with_breit=with_breit, npes=npes,
                        filename=new, mode=mode, IsComplex=True,
                    )
                    t_new = time.time() - t
                    self.assertLess(t_new, t_old, "%s breit=%d" % (mode, with_breit))
                    report.append(
                        "%s breit=%d: old %.2fs new %.2fs (%.2fx)"
                        % (mode, with_breit, t_old, t_new, t_old / t_new)
                    )
            print("\n".join(report))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    unittest.main(verbosity=2)
