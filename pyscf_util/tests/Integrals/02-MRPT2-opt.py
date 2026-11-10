#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for the optimized MRPT2 FCIDUMP writers and the binary format.

Run with::

    python pyscf_util/tests/Integrals/02-MRPT2-opt.py

Set ``RUN_SLOW=1`` to additionally run the H2O/cc-pVQZ performance checks.

The old writers in ``integral_MRPT2.py`` / ``integral_MRPT2_incore_fast.py``
are used as the oracle.  They are slow on purpose, but all systems below are
small (<= 28 orbitals) except the optional slow benchmark.
"""

import os
import re
import shutil
import tempfile
import time
import unittest
from collections import defaultdict

import numpy as np
import pyscf
from pyscf import ao2mo, mcscf, tools

from pyscf_util.Integrals import integral_MRPT2 as old_mod
from pyscf_util.Integrals import integral_MRPT2_incore_fast as old_fast_mod
from pyscf_util.Integrals import integral_MRPT2_binary as bin_mod
from pyscf_util.Integrals import integral_MRPT2_opt as opt_mod
from pyscf_util.misc.misc import _combine2, _combine4

RUN_SLOW = os.environ.get("RUN_SLOW", "0") not in ("0", "", "False", "false")

WATER = "O 0 0 0; H 0 0.757 0.586; H 0 -0.757 0.586"


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _pair(a, b):
    return (a, b) if a >= b else (b, a)


def _key8(a, b, c, d):
    """Pair based 8-fold key (1- or 0-based, order does not matter)."""
    ij = _pair(a, b)
    kl = _pair(c, d)
    return (ij, kl) if ij >= kl else (kl, ij)


def parse_fcidump(filename):
    """Return ``(records, h1e, ecore, orbsym, nelec, norb)``.

    ``records`` maps the 8-fold key to a list of values (duplicates kept so
    that the old writer's redundant records can be detected).
    """
    records = defaultdict(list)
    h1e = {}
    ecore = None
    orbsym = None
    nelec = None
    norb = None
    with open(filename) as fin:
        lines = fin.read().splitlines()
    end = 0
    for end, line in enumerate(lines):
        if "&END" in line:
            break
    head = "\n".join(lines[: end + 1])
    match = re.search(r"NORB\s*=\s*(\d+)", head)
    if match:
        norb = int(match.group(1))
    match = re.search(r"NELEC\s*=\s*(\d+)", head)
    if match:
        nelec = int(match.group(1))
    match = re.search(r"ORBSYM\s*=\s*([0-9,\s]+)", head)
    if match:
        orbsym = [int(x) for x in re.split(r"[,\s]+", match.group(1).strip()) if x]
    for line in lines[end + 1 :]:
        fields = line.split()
        if len(fields) < 5:
            continue
        value = float(fields[0])
        i, j, k, l = (int(x) for x in fields[1:5])
        if k == 0:
            if j == 0:
                ecore = value
            else:
                h1e[(i, j)] = value
            continue
        records[_key8(i, j, k, l)].append(value)
    return records, h1e, ecore, orbsym, nelec, norb


def binary_key8(records):
    """Set of 8-fold keys of a structured binary record array (1-based).

    The binary format stores 0-based indices; the keys are shifted back to the
    FCIDUMP convention so that they can be compared with ``parse_fcidump``.
    """
    keys = set()
    for p, q, r, s in zip(
        records["p"].tolist(),
        records["q"].tolist(),
        records["r"].tolist(),
        records["s"].tolist(),
    ):
        keys.add(_key8(p + 1, q + 1, r + 1, s + 1))
    return keys


def read_head_no_convert(filename):
    return tools.fcidump.read(filename, molpro_orbsym=False)


class _MoleculeCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp(prefix="mrpt2_opt_test_")
        cls.cases = {}
        for name, atom, basis, symmetry, split in (
            ("h2o_631g", WATER, "6-31g", "C2v", (1, 4)),
            ("n2_ccpvdz", "N 0 0 0; N 0 0 1.1", "cc-pvdz", "D2h", (2, 6)),
        ):
            mol = pyscf.gto.Mole()
            mol.atom = atom
            mol.basis = basis
            mol.symmetry = symmetry
            mol.verbose = 0
            mol.build()
            mf = pyscf.scf.RHF(mol).run()
            nfzc, nact = split
            cls.cases[name] = (mol, mf, mf.mo_coeff, nfzc, nact, mol.nao - nfzc - nact)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def path(self, name):
        return os.path.join(self.tmpdir, name)


# ---------------------------------------------------------------------------
# 1. ASCII writers: equivalence with the old implementations
# ---------------------------------------------------------------------------


class TestAsciiWriters(_MoleculeCase):
    def test_incore_byte_identical(self):
        for name, (mol, mf, mo, nfzc, nact, nvir) in self.cases.items():
            ref = self.path(name + "_incore_ref.txt")
            new = self.path(name + "_incore_new.txt")
            old_mod.fcidump_mrpt2(mol, mf, mo, nfzc, nact, nvir, ref, 1e-10)
            stats = opt_mod.fcidump_mrpt2_opt(
                mol, mf, mo, nfzc, nact, nvir, new, 1e-10
            )
            with open(ref, "rb") as f1, open(new, "rb") as f2:
                self.assertEqual(f1.read(), f2.read(), name)
            self.assertGreater(stats["nrecords"], 0)

    def test_outcore_byte_identical(self):
        for name, (mol, mf, mo, nfzc, nact, nvir) in self.cases.items():
            ref = self.path(name + "_outcore_ref.txt")
            new = self.path(name + "_outcore_new.txt")
            old_mod.fcidump_mrpt2_outcore(mol, mf, mo, nfzc, nact, nvir, ref, 1e-10)
            opt_mod.fcidump_mrpt2_outcore_opt(
                mol, mf, mo, nfzc, nact, nvir, new, 1e-10
            )
            with open(ref, "rb") as f1, open(new, "rb") as f2:
                self.assertEqual(f1.read(), f2.read(), name)

    def test_outcore_reuses_erifile(self):
        name = "h2o_631g"
        mol, mf, mo, nfzc, nact, nvir = self.cases[name]
        h5file = self.path("reuse_eri.h5")
        first = self.path("reuse_a.txt")
        second = self.path("reuse_b.txt")
        opt_mod.fcidump_mrpt2_outcore_opt(
            mol, mf, mo, nfzc, nact, nvir, first, 1e-10, erifile=h5file
        )
        self.assertTrue(os.path.exists(h5file))
        stats = opt_mod.fcidump_mrpt2_outcore_opt(
            mol, mf, mo, nfzc, nact, nvir, second, 1e-10, erifile=h5file
        )
        with open(first, "rb") as f1, open(second, "rb") as f2:
            self.assertEqual(f1.read(), f2.read())
        self.assertGreater(stats["nrecords"], 0)

    def test_fast_key_set_and_values(self):
        for name, (mol, mf, mo, nfzc, nact, nvir) in self.cases.items():
            for full in (False, True):
                ref = self.path("%s_fast_ref_%d.txt" % (name, full))
                new = self.path("%s_fast_new_%d.txt" % (name, full))
                old_fast_mod.fcidump_mrpt2_incore_fast(
                    mol, mf, mo, nfzc, nact, nvir, ref, 1e-10, full_integrals=full
                )
                opt_mod.fcidump_mrpt2_incore_fast_opt(
                    mol, mf, mo, nfzc, nact, nvir, new, 1e-10, full_integrals=full
                )
                rec_ref, _, _, _, _, _ = parse_fcidump(ref)
                rec_new, _, _, _, _, _ = parse_fcidump(new)
                self.assertEqual(
                    set(rec_ref), set(rec_new), "%s full=%s" % (name, full)
                )
                self.assertEqual(
                    sum(len(v) for v in rec_new.values()),
                    len(rec_new),
                    "duplicated records in the new fast writer",
                )
                for key, value in rec_new.items():
                    self.assertLessEqual(
                        abs(value[0] - rec_ref[key][0]), 1e-10, (name, key)
                    )
                for pair_a, pair_b in rec_new:
                    # parsed keys are 1-based
                    for idx in pair_a + pair_b:
                        self.assertLessEqual(idx, nfzc + nact + nvir)

    def test_hcore_and_ecore_identical(self):
        name = "h2o_631g"
        mol, mf, mo, nfzc, nact, nvir = self.cases[name]
        ref = self.path("hcore_ref.txt")
        new = self.path("hcore_new.txt")
        old_mod.fcidump_mrpt2(mol, mf, mo, nfzc, nact, nvir, ref, 1e-10)
        opt_mod.fcidump_mrpt2_opt(mol, mf, mo, nfzc, nact, nvir, new, 1e-10)
        _, h1_ref, e_ref, sym_ref, n_ref, _ = parse_fcidump(ref)
        _, h1_new, e_new, sym_new, n_new, _ = parse_fcidump(new)
        self.assertEqual(h1_ref, h1_new)
        self.assertAlmostEqual(e_ref, e_new, places=12)
        self.assertEqual(sym_ref, sym_new)
        self.assertEqual(n_ref, n_new)

    def test_matches_pyscf_reference(self):
        """The sparse writers must agree with pyscf's own FCIDUMP writer."""
        name = "h2o_631g"
        mol, mf, mo, nfzc, nact, nvir = self.cases[name]
        new = self.path("ref_cmp_new.txt")
        opt_mod.fcidump_mrpt2_opt(mol, mf, mo, nfzc, nact, nvir, new, 1e-10)
        rec_new, h1e_new, _, _, _, _ = parse_fcidump(new)
        ref = self.path("ref_cmp_pyscf.txt")
        tools.fcidump.from_mo(mol, ref, mo, tol=1e-10)
        rec_ref, h1e_ref, _, _, _, _ = parse_fcidump(ref)
        self.assertTrue(set(rec_new) <= set(rec_ref))
        for key, value in rec_new.items():
            self.assertLessEqual(abs(value[0] - rec_ref[key][0]), 1e-9, key)
        for (i, j), value in h1e_new.items():
            self.assertLessEqual(abs(value - h1e_ref[(i, j)]), 1e-12)


# ---------------------------------------------------------------------------
# 2. binary format
# ---------------------------------------------------------------------------


class TestBinaryFormat(_MoleculeCase):
    def _ascii_and_binary(self, name, mode="mrpt2"):
        mol, mf, mo, nfzc, nact, nvir = self.cases[name]
        ascii_file = self.path("%s_bin_%s.txt" % (name, mode))
        bin_file = self.path("%s_%s.bin" % (name, mode))
        if mode == "full":
            bin_mod.dump_mrpt2_binary(
                mol, mf, mo, nfzc, nact, nvir, bin_file, 1e-10, mode="full"
            )
            ref = self.path("%s_full_ref.txt" % name)
            tools.fcidump.from_mo(mol, ref, mo, tol=1e-10)
            ascii_file = ref
        else:
            old_mod.fcidump_mrpt2(mol, mf, mo, nfzc, nact, nvir, ascii_file, 1e-10)
            bin_mod.dump_mrpt2_binary(
                mol, mf, mo, nfzc, nact, nvir, bin_file, 1e-10, mode="mrpt2"
            )
        return ascii_file, bin_file

    def test_header_and_roundtrip(self):
        for name in self.cases:
            ascii_file, bin_file = self._ascii_and_binary(name)
            rec_ascii, h1e_ascii, ecore_ascii, _, nelec_ascii, norb_ascii = (
                parse_fcidump(ascii_file)
            )
            header = bin_mod.read_mrpt2_binary_header(bin_file)
            data = bin_mod.load_mrpt2_binary(bin_file)

            self.assertEqual(header["fmt"], bin_mod.MRPT2BIN_FMT)
            self.assertEqual(header["idxbase"], 0)
            self.assertEqual(header["nmo"], norb_ascii)
            self.assertEqual(data["nint"], header["nint"])
            self.assertEqual(header["nint"], len(data["eri"]))
            self.assertEqual(header["nelec"], nelec_ascii)
            self.assertAlmostEqual(data["ecore"], ecore_ascii, places=12)

            # indices must be identical (0-based in binary, 1-based in ASCII)
            self.assertEqual(binary_key8(data["eri"]), set(rec_ascii))
            for record in data["eri"]:
                key = _key8(
                    int(record["p"]) + 1,
                    int(record["q"]) + 1,
                    int(record["r"]) + 1,
                    int(record["s"]) + 1,
                )
                self.assertLessEqual(abs(record["v"] - rec_ascii[key][0]), 1e-12)
            # canonical form
            for record in data["eri"]:
                p, q, r, s = (int(record[x]) for x in "pqrs")
                self.assertGreaterEqual(p, q)
                self.assertGreaterEqual(r, s)
                self.assertGreaterEqual(
                    _combine2(p, q), _combine2(r, s), msg="not canonical"
                )
            # h1e
            for (i, j), value in h1e_ascii.items():
                self.assertLessEqual(abs(data["h1e"][i - 1, j - 1] - value), 1e-12)
            # no duplicate 8-fold keys
            self.assertEqual(len(set(binary_key8(data["eri"]))), len(data["eri"]))
            # size advantage
            self.assertLess(os.path.getsize(bin_file), 0.55 * os.path.getsize(ascii_file))

    def test_outcore_binary_matches_incore(self):
        mol, mf, mo, nfzc, nact, nvir = self.cases["h2o_631g"]
        f_inc = self.path("bin_inc.bin")
        f_out = self.path("bin_out.bin")
        bin_mod.dump_mrpt2_binary(
            mol, mf, mo, nfzc, nact, nvir, f_inc, 1e-10, mode="mrpt2"
        )
        bin_mod.dump_mrpt2_binary_outcore(
            mol, mf, mo, nfzc, nact, nvir, f_out, 1e-10, mode="mrpt2"
        )
        with open(f_inc, "rb") as f1, open(f_out, "rb") as f2:
            self.assertEqual(f1.read(), f2.read())

    def test_full_mode_against_pyscf(self):
        name = "h2o_631g"
        ascii_file, bin_file = self._ascii_and_binary(name, mode="full")
        converted = self.path("full_converted.txt")
        bin_mod.mrpt2_binary_to_fcidump(bin_file, converted)
        ref = read_head_no_convert(ascii_file)
        got = read_head_no_convert(converted)
        self.assertLessEqual(np.abs(ref["H1"] - got["H1"]).max(), 1e-12)
        self.assertLessEqual(np.abs(ref["H2"] - got["H2"]).max(), 1e-12)
        self.assertAlmostEqual(ref["ECORE"], got["ECORE"], places=12)
        self.assertEqual(ref["NORB"], got["NORB"])
        self.assertEqual(ref["NELEC"], got["NELEC"])

    def test_pyscf_dict_matches_fcidump_read(self):
        _, bin_file = self._ascii_and_binary("h2o_631g")
        converted = self.path("dict_converted.txt")
        bin_mod.mrpt2_binary_to_fcidump(bin_file, converted)
        ref = read_head_no_convert(converted)
        got = bin_mod.mrpt2_binary_to_pyscf_dict(bin_file)
        self.assertLessEqual(np.abs(ref["H1"] - got["H1"]).max(), 1e-12)
        self.assertLessEqual(np.abs(ref["H2"] - got["H2"]).max(), 1e-12)
        self.assertAlmostEqual(ref["ECORE"], got["ECORE"], places=12)
        for key in ("NORB", "NELEC", "MS2", "ORBSYM", "ISYM"):
            self.assertEqual(ref[key], got[key])

    def test_fast_binary_matches_old_fast(self):
        for name in self.cases:
            mol, mf, mo, nfzc, nact, nvir = self.cases[name]
            ref = self.path(name + "_fastbin_ref.txt")
            old_fast_mod.fcidump_mrpt2_incore_fast(
                mol, mf, mo, nfzc, nact, nvir, ref, 1e-10
            )
            rec_ref, _, _, _, _, _ = parse_fcidump(ref)
            bin_file = self.path(name + "_fast.bin")
            bin_mod.dump_mrpt2_binary_fast(
                mol, mf, mo, nfzc, nact, nvir, bin_file, 1e-10
            )
            data = bin_mod.load_mrpt2_binary(bin_file)
            self.assertEqual(binary_key8(data["eri"]), set(rec_ref))
            for record in data["eri"]:
                key = _key8(
                    int(record["p"]) + 1,
                    int(record["q"]) + 1,
                    int(record["r"]) + 1,
                    int(record["s"]) + 1,
                )
                self.assertLessEqual(abs(record["v"] - rec_ref[key][0]), 1e-10)

    def test_iter_mrpt2_binary(self):
        _, bin_file = self._ascii_and_binary("h2o_631g")
        chunks = list(bin_mod.iter_mrpt2_binary(bin_file, chunk=100))
        self.assertGreater(len(chunks), 1)
        total = sum(len(c) for c in chunks)
        self.assertEqual(total, bin_mod.read_mrpt2_binary_header(bin_file)["nint"])


# ---------------------------------------------------------------------------
# 3. edge cases
# ---------------------------------------------------------------------------


class TestEdgeCases(_MoleculeCase):
    def test_no_symmetry_molecule(self):
        mol = pyscf.gto.Mole()
        mol.atom = "O 0 0 0; H 0.757 0 0.586; H -0.757 0 0.586"
        mol.basis = "6-31g"
        mol.verbose = 0
        mol.build()
        self.assertIsNone(mol.irrep_name)
        mf = pyscf.scf.RHF(mol).run()
        mo = mf.mo_coeff
        nfzc, nact = 1, 4
        nvir = mol.nao - nfzc - nact
        with self.assertRaises(TypeError):
            old_mod.fcidump_mrpt2(mol, mf, mo, nfzc, nact, nvir, self.path("c1_o.txt"))
        out = self.path("c1_n.txt")
        opt_mod.fcidump_mrpt2_opt(mol, mf, mo, nfzc, nact, nvir, out, 1e-10)
        _, _, _, orbsym, _, _ = parse_fcidump(out)
        self.assertEqual(orbsym, [1] * mol.nao)

    def test_nao_larger_than_nmo(self):
        mol, mf, mo, nfzc, nact, _ = self.cases["h2o_631g"]
        nvir = 3  # nmo < nao: mo_coeff still has the full nao columns
        nmo = nfzc + nact + nvir
        self.assertLess(nmo, mol.nao)
        with self.assertRaises(Exception):
            old_fast_mod.fcidump_mrpt2_incore_fast(
                mol, mf, mo, nfzc, nact, nvir, self.path("small_o.txt")
            )
        out = self.path("small_n.txt")
        stats = opt_mod.fcidump_mrpt2_incore_fast_opt(
            mol, mf, mo, nfzc, nact, nvir, out, 1e-10
        )
        rec, _, _, _, _, _ = parse_fcidump(out)
        self.assertGreater(stats["nrecords"], 0)
        for pair_a, pair_b in rec:
            for idx in pair_a + pair_b:
                self.assertLessEqual(idx, nmo)

    def test_nfzc_zero_and_nvir_zero(self):
        mol, mf, mo, _, nact, _ = self.cases["h2o_631g"]
        nmo = mol.nao
        for nfzc, nvir in ((0, nmo - nact), (1, 0)):
            out = self.path("edge_%d_%d.txt" % (nfzc, nvir))
            stats = opt_mod.fcidump_mrpt2_opt(
                mol, mf, mo, nfzc, nact, nvir, out, 1e-10
            )
            self.assertGreaterEqual(stats["nrecords"], 0)
            rec, _, ecore, _, _, _ = parse_fcidump(out)
            self.assertIsNotNone(ecore)
            # binary must work as well
            bin_file = self.path("edge_%d_%d.bin" % (nfzc, nvir))
            bin_mod.dump_mrpt2_binary(
                mol, mf, mo, nfzc, nact, nvir, bin_file, 1e-10
            )
            data = bin_mod.load_mrpt2_binary(bin_file)
            self.assertEqual(len(data["eri"]), len(rec))

    def test_empty_record_set(self):
        mol, mf, mo, nfzc, nact, nvir = self.cases["h2o_631g"]
        out = self.path("empty.txt")
        stats = opt_mod.fcidump_mrpt2_opt(
            mol, mf, mo, nfzc, nact, nvir, out, 1e6
        )
        self.assertEqual(stats["nrecords"], 0)
        rec, _, ecore, _, _, _ = parse_fcidump(out)
        self.assertEqual(len(rec), 0)
        self.assertIsNotNone(ecore)
        bin_file = self.path("empty.bin")
        bin_mod.dump_mrpt2_binary(mol, mf, mo, nfzc, nact, nvir, bin_file, 1e6)
        data = bin_mod.load_mrpt2_binary(bin_file)
        self.assertEqual(data["nint"], 0)
        self.assertEqual(len(data["eri"]), 0)

    def test_uint16_limit(self):
        with self.assertRaises(ValueError):
            bin_mod.write_mrpt2_binary(
                self.path("toolarge.bin"),
                np.zeros((2, 2)),
                [],
                nmo=bin_mod.MRPT2_MAX_NMO + 1,
                nelec=2,
            )

    def test_truncated_file(self):
        mol, mf, mo, nfzc, nact, nvir = self.cases["h2o_631g"]
        good = self.path("trunc_good.bin")
        bin_mod.dump_mrpt2_binary(mol, mf, mo, nfzc, nact, nvir, good, 1e-10)
        bad = self.path("trunc_bad.bin")
        with open(good, "rb") as fin, open(bad, "wb") as fout:
            fout.write(fin.read()[:-8])
        with self.assertRaises(ValueError):
            bin_mod.read_mrpt2_binary_header(bad)

    def test_bad_magic(self):
        mol, mf, mo, nfzc, nact, nvir = self.cases["h2o_631g"]
        good = self.path("magic_good.bin")
        bin_mod.dump_mrpt2_binary(mol, mf, mo, nfzc, nact, nvir, good, 1e-10)
        header = bin_mod.read_mrpt2_binary_header(good)
        bad = self.path("magic_bad.bin")
        with open(good, "rb") as fin:
            payload = bytearray(fin.read())
        payload[header["head_end"] : header["head_end"] + 8] = b"XXXXXXXX"
        with open(bad, "wb") as fout:
            fout.write(payload)
        with self.assertRaises(ValueError):
            bin_mod.read_mrpt2_binary_header(bad)


# ---------------------------------------------------------------------------
# 4. the symmetry / index conventions the optimizations rely on
# ---------------------------------------------------------------------------


class TestSymmetryAssumptions(_MoleculeCase):
    def test_index_conventions(self):
        rng = np.random.default_rng(7)
        for _ in range(200):
            a, b, c, d = (int(x) for x in rng.integers(0, 40, 4))
            self.assertEqual(int(opt_mod.combine2(a, b)), int(_combine2(a, b)))
            self.assertEqual(
                int(opt_mod.combine4(a, b, c, d)), int(_combine4(a, b, c, d))
            )

    def test_canonicalize4(self):
        rng = np.random.default_rng(11)
        for _ in range(200):
            i, j, k, l = (int(x) for x in rng.integers(0, 30, 4))
            p, q, r, s = (
                int(np.asarray(x).ravel()[0])
                for x in opt_mod.canonicalize4(
                    np.array([i]), np.array([j]), np.array([k]), np.array([l])
                )
            )
            self.assertGreaterEqual(p, q)
            self.assertGreaterEqual(r, s)
            self.assertGreaterEqual(_combine2(p, q), _combine2(r, s))
            self.assertEqual(
                int(_combine4(i, j, k, l)), int(_combine4(p, q, r, s))
            )

    def test_aosym_s4_equals_s1(self):
        """4-fold AO symmetry must not change any block value."""
        mol, mf, mo, nfzc, nact, _ = self.cases["h2o_631g"]
        mc = mo[:, :nfzc]
        ma = mo[:, nfzc : nfzc + nact]
        mv = mo[:, nfzc + nact :]
        blocks = (
            (ma, ma, ma, ma),
            (ma, ma, mo, ma),
            (mv, ma, mv, ma),
            (ma, mc, ma, mc),
            (ma, ma, mc, mv),
            (ma, mc, ma, mv),
        )
        for block in blocks:
            ref = np.asarray(ao2mo.general(mol, block, aosym="s1", compact=False))
            sub = np.asarray(ao2mo.general(mol, block, aosym="s4", compact=False))
            self.assertEqual(ref.shape, sub.shape)
            self.assertLessEqual(np.abs(ref - sub).max(), 1e-12)

    def test_merged_pprs_prps_against_reference(self):
        mol, mf, mo, _, _, _ = self.cases["h2o_631g"]
        nao = mol.nao
        full = ao2mo.restore(1, np.asarray(ao2mo.full(mol, mo, aosym="4")).copy(), nao)
        pprs, prps = opt_mod._generate_eri_pprs_prps(mol, mo, verbose=0)
        self.assertEqual(pprs.shape, (nao, nao, nao))
        self.assertEqual(prps.shape, (nao, nao, nao))
        rng = np.random.default_rng(3)
        for _ in range(300):
            p, r, s = (int(x) for x in rng.integers(0, nao, 3))
            self.assertLessEqual(abs(pprs[p, r, s] - full[p, p, r, s]), 1e-10)
            self.assertLessEqual(abs(prps[p, r, s] - full[p, r, p, s]), 1e-10)
        # both target tensors are symmetric in r <-> s, which is what the
        # k <= l mirroring of the generator relies on.
        self.assertLessEqual(np.abs(pprs - pprs.transpose(0, 2, 1)).max(), 1e-10)
        self.assertLessEqual(np.abs(prps - prps.transpose(0, 2, 1)).max(), 1e-10)


# ---------------------------------------------------------------------------
# 5. generalized Fock
# ---------------------------------------------------------------------------


class TestGeneralizedFock(_MoleculeCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        mol, mf, mo, _, _, _ = cls.cases["h2o_631g"]
        cls.mc = mcscf.CASSCF(mf, 4, 4).run()
        cls.rdm1 = cls.mc.fcisolver.make_rdm1(
            cls.mc.ci, cls.mc.ncas, cls.mc.nelecas
        )
        cls.mo = cls.mc.mo_coeff

    def test_default_path(self):
        for kwargs in ({}, {"with_g3": True, "only_g3": True}, {"with_g3": True}):
            ref = old_mod.get_generalized_fock(self.mc, self.mo, self.rdm1, **kwargs)
            got = opt_mod.get_generalized_fock_opt(
                self.mc, self.mo, self.rdm1, **kwargs
            )
            self.assertLessEqual(np.abs(ref - got).max(), 1e-12, kwargs)

    def test_fallback_path(self):
        ref = old_mod.get_generalized_fock(
            self.mc, self.mo, self.rdm1, True, True, False
        )
        got = opt_mod.get_generalized_fock_opt(
            self.mc, self.mo, self.rdm1, True, True, False
        )
        self.assertLessEqual(np.abs(ref - got).max(), 1e-12)


# ---------------------------------------------------------------------------
# 6. optional performance check
# ---------------------------------------------------------------------------


@unittest.skipUnless(RUN_SLOW, "set RUN_SLOW=1 to run the cc-pVQZ benchmark")
class TestPerformance(unittest.TestCase):
    def test_large_case(self):
        mol = pyscf.gto.Mole()
        mol.atom = WATER
        mol.basis = "cc-pVQZ"
        mol.symmetry = "C2v"
        mol.verbose = 0
        mol.build()
        mf = pyscf.scf.RHF(mol).run()
        mo = mf.mo_coeff
        nfzc, nact = 1, 8
        nvir = mol.nao - nfzc - nact
        tmp = tempfile.mkdtemp(prefix="mrpt2_opt_perf_")
        try:
            report = []
            for tag, old_fn, new_fn in (
                ("incore", old_mod.fcidump_mrpt2, opt_mod.fcidump_mrpt2_opt),
                (
                    "outcore",
                    old_mod.fcidump_mrpt2_outcore,
                    opt_mod.fcidump_mrpt2_outcore_opt,
                ),
                (
                    "fast",
                    old_fast_mod.fcidump_mrpt2_incore_fast,
                    opt_mod.fcidump_mrpt2_incore_fast_opt,
                ),
            ):
                ref = os.path.join(tmp, tag + "_ref.txt")
                new = os.path.join(tmp, tag + "_new.txt")
                t0 = time.time()
                old_fn(mol, mf, mo, nfzc, nact, nvir, ref, 1e-10)
                t_old = time.time() - t0
                t0 = time.time()
                new_fn(mol, mf, mo, nfzc, nact, nvir, new, 1e-10)
                t_new = time.time() - t0
                report.append(
                    "%s: old %.2fs (%d bytes) new %.2fs (%d bytes)"
                    % (tag, t_old, os.path.getsize(ref), t_new, os.path.getsize(new))
                )
                self.assertLess(t_new, t_old, tag)
            print("\n".join(report))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    unittest.main(verbosity=2)
