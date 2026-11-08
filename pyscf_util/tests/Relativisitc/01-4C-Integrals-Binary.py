#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for the binary four-component integral dump.

Run with::

    python pyscf_util/tests/Relativisitc/01-4C-Integrals-Binary.py

What is checked
---------------
1. the binary record set is *identical* to the classic ASCII dump
   (in-core / out-of-core, real / complex, ``s4`` / ``s1``, Paldus/DG filters,
   ``no_2e``);
2. the ASCII file written next to the binary one keeps the usual FCIDUMP head
   (and nothing else by default) and its ``NORB`` / ``NELEC`` / ``ORBSYM``
   agree with the ASCII reference;
3. the byte layout is the one the C++ reader in
   ``iCIPT2_Loong/src/quantum_chemistry/integral_manager/Impl/integral_manager_ABC.h``
   expects: ``f8 value [, f8 imag] | u2 p | u2 q | u2 r | u2 s``, 1-based
   indices, 16/24 bytes per record, file size a multiple of the record size;
4. a small re-implementation of the C++ reader rebuilds the same ECORE,
   one- and two-electron integrals from the binary file;
5. ``binary_to_fcidump`` round-trips back to a readable ASCII FCIDUMP;
6. the D2H/Jz-adapted writer (:func:`FCIDUMP_Rela4C_SU2_binary`) -- the one
   feeding ``iCIPT2_Breit_D2H_ERICompressed`` -- agrees with the ASCII
   ``FCIDUMP_Rela4C_SU2_opt`` reference and stores real 16-byte records.
"""

import copy
import os
import re
import shutil
import struct
import tempfile
import unittest

import numpy as np

from pyscf import gto, scf

from pyscf_util.Relativisitc import integral_4C_binary as bin_mod
from pyscf_util.Relativisitc import integral_4C_opt as opt_mod

WATER = "O 0 0 0; H 0 -0.757 0.587; H 0 0.757 0.587"

HEAD_KEYS = ("NORB", "NELEC", "ORBSYM", "ISYM")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def parse_head(path):
    """Mimic the C++ head reader: NORB / NELEC / ORBSYM / ISYM + '&END'."""
    with open(path) as fin:
        lines = fin.readlines()
    head = {}
    n_integral_lines = 0
    ended = False
    for line in lines:
        if "&END" in line:
            ended = True
            continue
        if ended:
            if line.strip():
                n_integral_lines += 1
            continue
        if "&FCI" in line:
            ints = [int(x) for x in re.findall(r"\d+", line)]
            head["NORB"], head["NELEC"], head["MS2"] = ints[0], ints[1], ints[2]
        if "ORBSYM" in line:
            head["ORBSYM"] = [int(x) for x in re.findall(r"\d+", line)]
        if "ISYM" in line:
            head["ISYM"] = int(re.findall(r"\d+", line)[0])
    assert ended, "%s: no '&END' in the head" % path
    return head, n_integral_lines


def cpp_style_read(path, is_complex=False):
    """Re-implementation of the C++ binary branch of ``_read_fcidump``.

    Returns ``(records, h1e, ecore)`` with 1-based keys, exactly as the C++
    code dispatches them (``r == s == 0`` -> 1e/ECORE, else 2e; records with
    ``|v| < 1e-12`` are dropped).
    """
    with open(path, "rb") as fin:
        raw = fin.read()
    block = 24 if is_complex else 16
    assert len(raw) % block == 0, "file size is not a multiple of the block size"
    dtype = bin_mod.BINARY_RECORD_COMPLEX if is_complex else bin_mod.BINARY_RECORD_REAL
    data = np.frombuffer(raw, dtype=dtype)
    if is_complex:
        values = data["vr"].astype(np.complex128) + 1j * data["vi"]
    else:
        values = data["v"].astype(np.complex128)

    records = {}
    h1e = {}
    ecore = None
    for p, q, r, s, value in zip(
        data["p"].tolist(),
        data["q"].tolist(),
        data["r"].tolist(),
        data["s"].tolist(),
        values.tolist(),
    ):
        if abs(value) < 1e-12:
            continue
        if r == 0 and s == 0:
            if p == 0 and q == 0:
                ecore = value
            else:
                h1e[(p, q)] = value
        else:
            records[(p, q, r, s)] = value
    return records, h1e, ecore


def read_struct_record(path, offset=0, is_complex=False):
    """Decode one record with ``struct`` (independent of the numpy reader)."""
    with open(path, "rb") as fin:
        fin.seek(offset)
        raw = fin.read(24)
    if is_complex:
        re_, im, p, q, r, s = struct.unpack_from("<ddHHHH", raw, 0)
        return complex(re_, im), p, q, r, s
    v, p, q, r, s = struct.unpack_from("<dHHHH", raw, 0)
    return v, p, q, r, s


def cpp_epsilon_filter(records):
    """Drop the records the C++ reader would drop (``|v| < 1e-12``)."""
    return {key: value for key, value in records.items() if abs(value) >= 1e-12}


# ---------------------------------------------------------------------------
# plain four-component (no D2H adaptation)
# ---------------------------------------------------------------------------


class _FourComponentCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp(prefix="rela4c_binary_test_")
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

    def dump_pair(self, tag, mode="incore", is_complex=False, **kwargs):
        """ASCII reference (``FCIDUMP_Rela4C_opt``) + binary dump."""
        ref = self.path("ref_%s" % tag)
        head = self.path("bin_%s" % tag)
        opt_mod.FCIDUMP_Rela4C_opt(
            self.mol,
            self.mf,
            with_breit=True,
            npes=self.npes,
            filename=ref,
            mode=mode,
            IsComplex=is_complex,
            **kwargs,
        )
        bin_mod.FCIDUMP_Rela4C_binary(
            self.mol,
            self.mf,
            with_breit=True,
            npes=self.npes,
            filename=head,
            mode=mode,
            IsComplex=is_complex,
            **kwargs,
        )
        return ref, head, head + ".bin"


class TestBinaryVsAscii(_FourComponentCase):
    def _compare(self, ref, binary, is_complex, label):
        result = bin_mod.compare_binary_with_fcidump(
            binary, ref, is_complex=is_complex, tol=1e-10
        )
        self.assertGreater(result["n_2e"], 0, label)
        return result

    def test_incore_real(self):
        ref, _head, binary = self.dump_pair("incore_real", is_complex=False)
        stats = self._compare(ref, binary, False, "incore real")
        summary = bin_mod.binary_summary(binary, is_complex=False)
        self.assertEqual(summary["record_size"], 16)
        self.assertEqual(summary["n_2e"], stats["n_2e"])
        self.assertEqual(summary["n_1e"], stats["n_1e"])
        self.assertEqual(summary["n_ecore"], 1)

    def test_incore_complex(self):
        ref, _head, binary = self.dump_pair("incore_complex", is_complex=True)
        stats = self._compare(ref, binary, True, "incore complex")
        summary = bin_mod.binary_summary(binary, is_complex=True)
        self.assertEqual(summary["record_size"], 24)
        self.assertEqual(summary["n_2e"], stats["n_2e"])

    def test_outcore_real(self):
        ref, _head, binary = self.dump_pair("outcore_real", mode="outcore")
        self._compare(ref, binary, False, "outcore real")

    def test_outcore_complex(self):
        ref, _head, binary = self.dump_pair(
            "outcore_complex", mode="outcore", is_complex=True
        )
        self._compare(ref, binary, True, "outcore complex")

    def test_incore_equals_outcore(self):
        _, _, bin_in = self.dump_pair("io_in", mode="incore")
        _, _, bin_out = self.dump_pair("io_out", mode="outcore")
        rec_in = bin_mod.read_binary_records(bin_in, is_complex=False)
        rec_out = bin_mod.read_binary_records(bin_out, is_complex=False)
        self.assertEqual(len(rec_in), len(rec_out))
        np.testing.assert_array_equal(rec_in["p"], rec_out["p"])
        np.testing.assert_array_equal(rec_in["q"], rec_out["q"])
        np.testing.assert_array_equal(rec_in["r"], rec_out["r"])
        np.testing.assert_array_equal(rec_in["s"], rec_out["s"])
        self.assertLessEqual(float(np.max(np.abs(rec_in["v"] - rec_out["v"]))), 1e-12)

    def test_debug_full_dump(self):
        ref, _head, binary = self.dump_pair("debug_s1", is_complex=False, debug=True)
        stats = self._compare(ref, binary, False, "debug s1")
        # s1 must be a superset of the s4 (4-fold representative) record set
        ref_s4, _h, bin_s4 = self.dump_pair("debug_s4")
        stats_s4 = self._compare(ref_s4, bin_s4, False, "s4")
        self.assertGreater(stats["n_2e"], stats_s4["n_2e"])

    def test_paldus_and_dg_filters(self):
        for tag, kwargs in (
            ("pal123", {"retain_paldus_123": False}),
            ("pal456", {"retain_paldus_456": False}),
            ("onlydg", {"only_dg": True}),
        ):
            ref, _head, binary = self.dump_pair(tag, is_complex=False, **kwargs)
            stats = self._compare(ref, binary, False, tag)
            self.assertGreater(stats["n_2e"], 0)

    def test_no_2e(self):
        head = self.path("bin_no2e")
        bin_mod.FCIDUMP_Rela4C_binary(
            self.mol,
            self.mf,
            with_breit=True,
            npes=self.npes,
            filename=head,
            IsComplex=False,
            no_2e=True,
        )
        summary = bin_mod.binary_summary(head + ".bin", is_complex=False)
        self.assertEqual(summary["n_2e"], 0)
        self.assertGreater(summary["n_1e"], 0)
        self.assertEqual(summary["n_ecore"], 1)


class TestHeadAndLayout(_FourComponentCase):
    def test_ascii_file_is_head_only(self):
        ref, head, binary = self.dump_pair("headonly", is_complex=False)
        head_parsed, n_integral_lines = parse_head(head)
        ref_parsed, _ = parse_head(ref)
        self.assertEqual(n_integral_lines, 0)
        for key in HEAD_KEYS:
            self.assertEqual(head_parsed[key], ref_parsed[key], key)
        self.assertEqual(head_parsed["NORB"], self.npes // 2)

    def test_head_means_the_binary_holds_everything(self):
        _ref, head, binary = self.dump_pair("head_everything", is_complex=False)
        head_parsed, _n = parse_head(head)
        norb = head_parsed["NORB"]
        records, h1e, ecore = cpp_style_read(binary)
        # ECORE
        self.assertIsNotNone(ecore)
        self.assertAlmostEqual(ecore, self.mol.get_enuc(), places=12)
        # one-electron: shape implied by the head, i.e. 2*norb orbitals
        self.assertTrue(all(1 <= p <= 2 * norb for p, _q in h1e))
        self.assertTrue(all(1 <= q <= 2 * norb for _p, q in h1e))
        # two-electron: Coulomb block + Breit block (offset by 2*norb)
        self.assertTrue(all(1 <= p <= 4 * norb for p, _q, _r, _s in records))
        self.assertTrue(any(p > 2 * norb for p, _q, _r, _s in records))

    def test_cpp_reader_reproduces_ascii(self):
        ref, _head, binary = self.dump_pair("cpp_read", is_complex=False)
        rec_ref, h1_ref, ecore_ref, _head_text = bin_mod.parse_fcidump_records(ref)
        records, h1e, ecore = cpp_style_read(binary)
        # Both the ASCII and the binary reader drop records with |v| < 1e-12.
        # For an IsComplex=False dump only the real part is stored, so records
        # whose real part is tiny but whose magnitude is not are dropped by
        # the reader on *both* paths -- apply the same filter to the reference.
        rec_ref = cpp_epsilon_filter(rec_ref)
        h1_ref = cpp_epsilon_filter(h1_ref)
        self.assertEqual(set(records), set(rec_ref))
        self.assertEqual(set(h1e), set(h1_ref))
        self.assertLessEqual(max(abs(records[k] - rec_ref[k]) for k in rec_ref), 1e-10)
        self.assertLessEqual(max(abs(h1e[k] - h1_ref[k]) for k in h1_ref), 1e-10)
        self.assertAlmostEqual(ecore, ecore_ref.real, places=10)

    def test_struct_level_layout(self):
        _ref, _head, binary = self.dump_pair("layout", is_complex=False)
        size = os.path.getsize(binary)
        self.assertEqual(size % 16, 0)
        records = bin_mod.read_binary_records(binary, is_complex=False)
        # last record: the ECORE line written by dump_ecore_binary
        value, p, q, r, s = read_struct_record(binary, size - 16)
        self.assertAlmostEqual(value, self.mol.get_enuc(), places=12)
        self.assertEqual((p, q, r, s), (0, 0, 0, 0))
        self.assertEqual(int(records["p"].max()) <= 2 * self.npes, True)
        self.assertEqual(int(records["s"].max()) <= 2 * self.npes, True)

    def test_struct_level_layout_complex(self):
        _ref, _head, binary = self.dump_pair("layout_c", is_complex=True)
        size = os.path.getsize(binary)
        self.assertEqual(size % 24, 0)
        value, p, q, r, s = read_struct_record(binary, size - 24, is_complex=True)
        self.assertAlmostEqual(value.real, self.mol.get_enuc(), places=12)
        self.assertAlmostEqual(value.imag, 0.0, places=12)
        self.assertEqual((p, q, r, s), (0, 0, 0, 0))

    def test_write_ascii_integrals_option(self):
        ref, head, binary = self.dump_pair("withascii", is_complex=False)
        head_full = self.path("full_ascii")
        bin_mod.FCIDUMP_Rela4C_binary(
            self.mol,
            self.mf,
            with_breit=True,
            npes=self.npes,
            filename=head_full,
            IsComplex=False,
            write_ascii_integrals=True,
        )
        rec_ref, h1_ref, ecore_ref, _ = bin_mod.parse_fcidump_records(ref)
        rec_full, h1_full, ecore_full, _ = bin_mod.parse_fcidump_records(head_full)
        self.assertEqual(set(rec_ref), set(rec_full))
        self.assertEqual(set(h1_ref), set(h1_full))
        self.assertAlmostEqual(ecore_ref.real, ecore_full.real, places=12)
        # the binary file must still be the same as without the ASCII records
        self.assertEqual(os.path.getsize(binary), os.path.getsize(head_full + ".bin"))

    def test_binary_to_fcidump_roundtrip(self):
        ref, head, binary = self.dump_pair("roundtrip", is_complex=False)
        out = self.path("roundtrip.txt")
        bin_mod.binary_to_fcidump(binary, out, head_filename=head)
        rec_ref, h1_ref, ecore_ref, head_ref = bin_mod.parse_fcidump_records(ref)
        rec_out, h1_out, ecore_out, head_out = bin_mod.parse_fcidump_records(out)
        self.assertEqual(set(rec_ref), set(rec_out))
        self.assertEqual(set(h1_ref), set(h1_out))
        self.assertLessEqual(max(abs(rec_ref[k] - rec_out[k]) for k in rec_ref), 1e-10)
        self.assertAlmostEqual(ecore_ref.real, ecore_out.real, places=12)
        # the head is copied verbatim
        self.assertEqual(head_ref.split(), head_out.split())

    def test_fcidump_to_binary(self):
        """A classic ASCII FCIDUMP can be converted without recomputation."""
        ref, _head, binary = self.dump_pair("to_bin_ref", is_complex=False)
        head = self.path("to_bin_head")
        converted = self.path("to_bin.bin")
        stats = bin_mod.fcidump_to_binary(ref, head, converted)
        self.assertFalse(stats["is_complex"])
        self.assertEqual(stats["nrecords"], os.path.getsize(converted) // 16)

        # the converted records must be the same, in the same order (the values
        # are only as precise as the "%18.12E" ASCII text)
        original = bin_mod.read_binary_records(binary, is_complex=False)
        again = bin_mod.read_binary_records(converted, is_complex=False)
        np.testing.assert_array_equal(original["p"], again["p"])
        np.testing.assert_array_equal(original["q"], again["q"])
        np.testing.assert_array_equal(original["r"], again["r"])
        np.testing.assert_array_equal(original["s"], again["s"])
        np.testing.assert_allclose(original["v"], again["v"], rtol=1e-11, atol=1e-11)

        # and the head-only file is a valid FCIDUMP head
        head_parsed, n_lines = parse_head(head)
        ref_parsed, _ = parse_head(ref)
        self.assertEqual(n_lines, 0)
        for key in HEAD_KEYS:
            self.assertEqual(head_parsed[key], ref_parsed[key], key)

    def test_fcidump_to_binary_complex(self):
        ref, _head, binary = self.dump_pair("to_bin_c_ref", is_complex=True)
        head = self.path("to_bin_c_head")
        converted = self.path("to_bin_c.bin")
        stats = bin_mod.fcidump_to_binary(ref, head, converted)
        self.assertTrue(stats["is_complex"])
        original = bin_mod.read_binary_records(binary, is_complex=True)
        again = bin_mod.read_binary_records(converted, is_complex=True)
        np.testing.assert_allclose(original["vr"], again["vr"], rtol=1e-11, atol=1e-11)
        np.testing.assert_allclose(original["vi"], again["vi"], rtol=1e-11, atol=1e-11)

    def test_wrong_record_size_is_reported(self):
        # a 16-byte (real) file cannot be a multiple of the 24-byte records
        tiny = self.path("tiny.bin")
        with bin_mod.BinaryIntegralWriter(tiny, is_complex=False) as writer:
            writer.add_records([1], [2], [3], [4], [0.5])
        self.assertEqual(os.path.getsize(tiny), 16)
        with self.assertRaises(ValueError):
            bin_mod.read_binary_records(tiny, is_complex=True)
        with self.assertRaises(ValueError):
            bin_mod.binary_summary(tiny, is_complex=True)

    def test_index_limit(self):
        with self.assertRaises(ValueError):
            bin_mod._check_index_limit(bin_mod.MAX_ORBITAL_INDEX // 2 + 1)


# ---------------------------------------------------------------------------
# D2H / Jz-adapted writer (the iCIPT2_Breit_D2H input)
# ---------------------------------------------------------------------------


class TestSU2Binary(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp(prefix="rela4c_su2_binary_test_")
        cls.mol = gto.M(
            atom="F 0 0 0",
            basis="sto-3g",
            verbose=0,
            charge=-1,
            spin=0,
            symmetry="d2h",
        )
        cls.mol.build()
        cls.mf = scf.dhf.RDHF(cls.mol)
        cls.mf.conv_tol = 1e-10
        cls.mf.kernel()
        cls.mf.with_breit = True
        cls.mo0 = np.array(cls.mf.mo_coeff)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def path(self, name):
        return os.path.join(self.tmpdir, name)

    def _fresh_mf(self):
        mf = copy.copy(self.mf)
        mf.mo_coeff = self.mo0.copy()
        return mf

    def test_su2_binary_matches_su2_ascii(self):
        ref = self.path("su2_ref")
        head = self.path("su2_bin")
        _c, _b, ref_parity, _mo_ref = opt_mod.FCIDUMP_Rela4C_SU2_opt(
            self.mol, self._fresh_mf(), True, filename=ref, mode="incore"
        )
        _c, _b, mo_parity, _mo = bin_mod.FCIDUMP_Rela4C_SU2_binary(
            self.mol, self._fresh_mf(), True, filename=head, mode="incore"
        )
        result = bin_mod.compare_binary_with_fcidump(
            head + ".bin", ref, is_complex=False, tol=1e-10, verbose=True
        )
        self.assertGreater(result["n_2e"], 0)
        self.assertEqual(list(mo_parity), list(ref_parity))

        summary = bin_mod.binary_summary(head + ".bin", is_complex=False)
        self.assertEqual(summary["record_size"], 16)  # real double manager
        self.assertEqual(summary["n_ecore"], 1)

        head_parsed, n_lines = parse_head(head)
        ref_parsed, _ = parse_head(ref)
        self.assertEqual(n_lines, 0)
        for key in HEAD_KEYS:
            self.assertEqual(head_parsed[key], ref_parsed[key], key)

    def test_su2_binary_rejects_too_many_spinors(self):
        with self.assertRaises(ValueError):
            bin_mod._check_index_limit(40000)


if __name__ == "__main__":
    unittest.main(verbosity=2)
