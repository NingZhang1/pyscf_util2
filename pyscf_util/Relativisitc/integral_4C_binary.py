#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Binary four-component integral dump (FCIDUMP head + binary integrals).

Motivation
----------
The four-component FCIDUMP files produced by
:mod:`pyscf_util.Relativisitc.integral_4C` (and the faster
:mod:`pyscf_util.Relativisitc.integral_4C_opt`) are dominated by the
two-electron records: one `` %18.12E %4d %4d %4d %4d`` line is ~33 bytes and
even a small basis easily reaches tens of millions of records.  Both the
writing (Python formatting) and the reading (C++ ``atof`` loop) become the
bottleneck.

This module keeps the **ASCII head exactly as it is** -- so ``NORB``,
``NELEC``, ``MS2``, ``ORBSYM``, ``ISYM`` remain human readable and every tool
that only needs the head keeps working -- and stores **all integrals** as raw
fixed-size binary records, which is the format consumed by the Loong iCIPT2
codes (e.g. ``src/app/iCIPT2_Breit_D2H_ERICompressed.cpp``).

File layout
-----------
Two files are produced, e.g. ``FCIDUMP_F`` (ASCII) and ``FCIDUMP_F.bin``::

    FCIDUMP_F                         # head only by default
     &FCI NORB=  40,NELEC=36,MS2=0,
      ORBSYM=0,0,0,...
      ISYM=1,
     &END

    FCIDUMP_F.bin                     # raw records, no head, no padding
      real    record : f8 value | u2 p | u2 q | u2 r | u2 s   (16 bytes)
      complex record : f8 re | f8 im | u2 p | u2 q | u2 r | u2 s (24 bytes)

Conventions (identical to the ASCII FCIDUMP integral lines)
-----------------------------------------------------------
* indices are **1-based**; ``0`` is the "absent index" placeholder:
  ``r = s = 0`` marks a one-electron record, ``p = q = r = s = 0`` marks the
  nuclear repulsion (``ECORE``);
* the Breit/Gaunt records use the index offset ``+npes`` exactly like the
  ASCII dump, so Coulomb and Breit live in the same stream;
* values are little-endian ``float64`` (the imaginary part is only present
  when ``IsComplex=True``);
* the record set (thresholding, Paldus/DG filters, 4-fold or full symmetry)
  is exactly the same as the ASCII writer -- the very same filters and loops
  of :mod:`pyscf_util.Relativisitc.integral_4C_opt` are reused -- so the two
  files describe the same integrals; the binary file keeps the full
  ``float64`` precision while the ASCII file is rounded to the
  ``%18.12E`` print precision;
* orbital indices are ``uint16``: ``2 * npes <= 65535`` is required (the
  Breit records are offset by ``npes``).

Matching C++ reader
-------------------
``QuantumChemistryLoong::Integrals::IntegralManager_ABC::_read_fcidump``
reads the head from the ASCII file (2nd command line argument) and, when a
binary file is given (3rd argument), parses *all* integrals from it::

    size_t blockSize = sizeof(double) + 4 * sizeof(uint16_t);   // 16
    if (!is_real) blockSize += sizeof(double);                  // 24
    // memcpy: value [, imag], p, q, r, s  (16 bit unsigned each)
    // p = q = r = s = 0            -> energy_core / repnuc
    // r = s = 0                    -> one-electron integral
    // |value| < 1e-12              -> record dropped

Consequences:

* the binary file must contain the **one-electron integrals and the nuclear
  repulsion as well** (they are not read from the ASCII file any more);
* ``IsComplex`` of the dump must match the scalar type of the C++ integral
  manager: ``Integrals4C_D2h_limited_memory`` (used by the ``*_D2H`` apps)
  is ``double``, i.e. **real 16-byte records** -> use ``IsComplex=False``,
  which is what :func:`FCIDUMP_Rela4C_SU2_binary` does.

Typical usage
-------------
>>> from pyscf_util.Relativisitc.integral_4C_binary import (
...     FCIDUMP_Rela4C_SU2_binary,
... )
>>> FCIDUMP_Rela4C_SU2_binary(mol, mf, True, filename="FCIDUMP_F", mode="outcore")
# -> FCIDUMP_F      : head only
# -> FCIDUMP_F.bin  : all integrals (Coulomb + Breit + h1e + ECORE)
#
# ./iCIPT2_Breit_D2H_ERICompressed  input  FCIDUMP_F  FCIDUMP_F.bin

Reading / checking
------------------
::

    iter_binary_records / read_binary_records   # structured numpy arrays
    binary_summary                              # record counts per type
    binary_to_fcidump                           # binary -> classic ASCII
    fcidump_to_binary                           # classic ASCII -> binary
    parse_fcidump_records                       # ASCII -> dict
    compare_binary_with_fcidump                 # record-by-record check

Set ``write_ascii_integrals=True`` to *additionally* write the complete ASCII
FCIDUMP (useful to cross-check the binary file, see the test suite).
"""

from __future__ import annotations

import os
import shutil
import tempfile
import time

import numpy as np
from pyscf import lib, tools

from pyscf_util.Relativisitc.integral_4C_opt import (
    _BREIT_TERMS,
    _FMT_ECORE_C,
    _FMT_ECORE_R,
    _GAUNT_TERMS,
    _build_breit_incore,
    _build_breit_outcore,
    _build_coulomb_incore,
    _build_coulomb_outcore,
    _dump_h1e_opt,
    _filter_mask,
    _mo_parts,
    _read_row,
    dump_2e_opt,
    dump_2e_outcore_opt,
)

__all__ = [
    "BINARY_RECORD_REAL",
    "BINARY_RECORD_COMPLEX",
    "MAX_ORBITAL_INDEX",
    "BinaryIntegralWriter",
    "dump_2e_binary",
    "dump_2e_outcore_binary",
    "dump_h1e_binary",
    "dump_ecore_binary",
    "FCIDUMP_Rela4C_binary",
    "FCIDUMP_Rela4C_SU2_binary",
    "read_binary_records",
    "iter_binary_records",
    "binary_summary",
    "binary_to_fcidump",
    "fcidump_to_binary",
    "parse_fcidump_records",
    "compare_binary_with_fcidump",
]

# ---------------------------------------------------------------------------
# record layout
# ---------------------------------------------------------------------------

#: 16-byte record of a real integral (value first, then the four indices).
BINARY_RECORD_REAL = np.dtype(
    [("v", "<f8"), ("p", "<u2"), ("q", "<u2"), ("r", "<u2"), ("s", "<u2")]
)

#: 24-byte record of a complex integral.
BINARY_RECORD_COMPLEX = np.dtype(
    [
        ("vr", "<f8"),
        ("vi", "<f8"),
        ("p", "<u2"),
        ("q", "<u2"),
        ("r", "<u2"),
        ("s", "<u2"),
    ]
)

assert BINARY_RECORD_REAL.itemsize == 16
assert BINARY_RECORD_COMPLEX.itemsize == 24

#: ``uint16`` orbital indices; ``0`` is reserved for "absent index".
MAX_ORBITAL_INDEX = 0xFFFF

_ASCII_FMT_REAL = "%18.12E %4d %4d %4d %4d\n"
_ASCII_FMT_COMPLEX = "%18.12E %18.12E %4d %4d %4d %4d\n"

#: threshold used by the C++ reader: smaller records are silently dropped.
CPP_ZERO_THRESHOLD = 1e-12

_DEFAULT_BUFFER_RECORDS = 1 << 20


def _record_dtype(is_complex):
    return BINARY_RECORD_COMPLEX if is_complex else BINARY_RECORD_REAL


def _log(verbose, msg):
    if verbose > 0:
        print(msg)


def _check_index_limit(npes):
    """The Breit records are offset by ``npes`` -> ``2 * npes`` must fit."""
    if 2 * int(npes) > MAX_ORBITAL_INDEX:
        raise ValueError(
            "npes=%d cannot be stored with uint16 orbital indices: the Breit "
            "records are offset by npes, so 2*npes=%d must be <= %d"
            % (int(npes), 2 * int(npes), MAX_ORBITAL_INDEX)
        )


# ---------------------------------------------------------------------------
# writer
# ---------------------------------------------------------------------------


class BinaryIntegralWriter:
    """Buffered writer of raw four-component integral records.

    Parameters
    ----------
    filename : str
        Binary output file (``.bin``).  It is truncated on creation.
    is_complex : bool
        ``False`` -> 16-byte records, ``True`` -> 24-byte records.
    buffer_records : int
        Number of records buffered before a write; purely a performance knob.

    Notes
    -----
    ``add_records`` expects 1-based indices (``0`` = absent index) and real or
    complex values.  When the writer is real, the imaginary part of complex
    input is dropped, exactly like the ``IsComplex=False`` ASCII dump.
    """

    def __init__(
        self, filename, is_complex=False, buffer_records=_DEFAULT_BUFFER_RECORDS
    ):
        self.filename = filename
        self.is_complex = bool(is_complex)
        self.dtype = _record_dtype(self.is_complex)
        self.record_size = int(self.dtype.itemsize)
        self.buffer_records = int(buffer_records)
        if self.buffer_records < 1:
            raise ValueError("buffer_records must be positive")
        self.nrecords = 0
        self._fh = open(filename, "wb")
        self._buf = []
        self._buffered = 0

    # -- public API ------------------------------------------------------

    def add_records(self, p, q, r, s, values):
        """Append ``len(values)`` records; returns the number of records."""
        values = np.asarray(values).reshape(-1)
        count = int(values.size)
        if count == 0:
            return 0

        p = self._index_array(p, count, "p")
        q = self._index_array(q, count, "q")
        r = self._index_array(r, count, "r")
        s = self._index_array(s, count, "s")

        record = np.empty(count, dtype=self.dtype)
        record["p"] = p
        record["q"] = q
        record["r"] = r
        record["s"] = s
        if self.is_complex:
            values = np.asarray(values, dtype=np.complex128)
            record["vr"] = values.real
            record["vi"] = values.imag
        else:
            if np.iscomplexobj(values):
                # mirrors the ASCII "IsComplex=False" writer
                values = values.real
            record["v"] = np.asarray(values, dtype=np.float64)

        self._buf.append(record)
        self._buffered += count
        self.nrecords += count
        if self._buffered >= self.buffer_records:
            self.flush()
        return count

    def flush(self):
        if not self._buf:
            return
        if len(self._buf) == 1:
            payload = self._buf[0].tobytes()
        else:
            payload = np.concatenate(self._buf).tobytes()
        self._fh.write(payload)
        self._buf = []
        self._buffered = 0

    def close(self):
        if self._fh is None:
            return
        self.flush()
        self._fh.close()
        self._fh = None

    # -- helpers ---------------------------------------------------------

    @staticmethod
    def _index_array(arr, count, name):
        arr = np.asarray(arr, dtype=np.int64).reshape(-1)
        if arr.size != count:
            raise ValueError(
                "index array '%s' has %d entries, expected %d" % (name, arr.size, count)
            )
        if count:
            low = int(arr.min())
            high = int(arr.max())
            if low < 0 or high > MAX_ORBITAL_INDEX:
                raise ValueError(
                    "orbital index '%s' out of uint16 range (1-based, 0 allowed): "
                    "%d .. %d" % (name, low, high)
                )
        return arr

    # -- context manager / bookkeeping -----------------------------------

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False

    @property
    def nbytes(self):
        return int(self.nrecords) * self.record_size

    def __repr__(self):
        return "BinaryIntegralWriter(%r, is_complex=%s, nrecords=%d)" % (
            self.filename,
            self.is_complex,
            self.nrecords,
        )


# ---------------------------------------------------------------------------
# 2e records -- identical filters/loops as the ASCII writers
# ---------------------------------------------------------------------------


def _emit_row_binary(
    writer,
    i,
    j,
    row_c,
    row_b,
    tol,
    retain_123,
    retain_456,
    only_dg,
    offset,
):
    """Emit the binary records of one ``(i, j)`` row.

    ``row_c`` / ``row_b`` have shape ``(nk, nl)`` where ``k`` covers the
    symmetry-allowed range and ``l`` all orbitals -- exactly the slices used
    by :func:`pyscf_util.Relativisitc.integral_4C_opt._emit_row`, so the record
    set *and the record order* are the same as the ASCII dump.
    """
    nk, nl = row_c.shape
    allowed = _filter_mask(
        i // 2,
        j // 2,
        (np.arange(nk) // 2)[:, None],
        (np.arange(nl) // 2)[None, :],
        retain_123,
        retain_456,
        only_dg,
    )
    mask_c = (np.abs(row_c) > tol) & allowed
    if row_b is not None:
        mask_b = (np.abs(row_b) > tol) & allowed
        union = mask_c | mask_b
    else:
        mask_b = None
        union = mask_c

    k_idx, l_idx = np.nonzero(union)
    if k_idx.size == 0:
        return

    k1 = k_idx + 1
    l1 = l_idx + 1
    i1 = np.full(k_idx.size, i + 1, dtype=np.int64)
    j1 = np.full(k_idx.size, j + 1, dtype=np.int64)

    if row_b is None:
        writer.add_records(i1, j1, k1, l1, row_c[k_idx, l_idx])
        return

    # interleave Coulomb and Breit records so that the order matches the
    # ASCII writer: c(k0,l0), b(k0,l0), c(k0,l1), b(k0,l1), ...
    keep_c = mask_c[k_idx, l_idx]
    keep_b = mask_b[k_idx, l_idx]
    n = int(k1.size)
    keep = np.empty(2 * n, dtype=bool)
    keep[0::2] = keep_c
    keep[1::2] = keep_b
    selection = np.nonzero(keep)[0]
    if selection.size == 0:
        return

    dtype = np.result_type(row_c.dtype, row_b.dtype)
    values = np.empty(2 * n, dtype=dtype)
    values[0::2] = row_c[k_idx, l_idx]
    values[1::2] = row_b[k_idx, l_idx]

    # the Breit/Gaunt records carry the "+npes" orbital offset on all four
    # indices, exactly like the ASCII writer
    shift = np.zeros(2 * n, dtype=np.int64)
    shift[1::2] = int(offset)
    pp = (np.repeat(i1, 2) + shift)[selection]
    qq = (np.repeat(j1, 2) + shift)[selection]
    rr = (np.repeat(k1, 2) + shift)[selection]
    ss = (np.repeat(l1, 2) + shift)[selection]

    writer.add_records(pp, qq, rr, ss, values[selection])


def dump_2e_binary(
    writer,
    int2e_coulomb,
    int2e_breit,
    with_breit,
    symmetry="s4",
    tol=1e-8,
    retain_paldus_123=True,
    retain_paldus_456=True,
    only_dg=False,
):
    """In-core binary 2e dump; same record set as ``dump_2e_opt``.

    Parameters
    ----------
    writer : BinaryIntegralWriter
    int2e_coulomb, int2e_breit : ndarray
        ``(npes, npes, npes, npes)`` MO integrals (``int2e_breit`` may be
        ``None`` when ``with_breit`` is false).
    symmetry : {"s4", "s1"}
        ``s4`` keeps ``p >= q`` and ``r <= p`` (4-fold representative set);
        ``s1`` writes every element (used for the ``debug`` files).
    """
    npes = int2e_coulomb.shape[0]
    _check_index_limit(npes)
    if symmetry == "s4":
        for i in range(npes):
            for j in range(i + 1):
                _emit_row_binary(
                    writer,
                    i,
                    j,
                    int2e_coulomb[i, j, : i + 1, :],
                    int2e_breit[i, j, : i + 1, :] if with_breit else None,
                    tol,
                    retain_paldus_123,
                    retain_paldus_456,
                    only_dg,
                    npes,
                )
    elif symmetry == "s1":
        for i in range(npes):
            for j in range(npes):
                _emit_row_binary(
                    writer,
                    i,
                    j,
                    int2e_coulomb[i, j, :, :],
                    int2e_breit[i, j, :, :] if with_breit else None,
                    tol,
                    retain_paldus_123,
                    retain_paldus_456,
                    only_dg,
                    npes,
                )
    else:
        raise ValueError("Unknown symmetry %s" % symmetry)
    return writer


def dump_2e_outcore_binary(
    writer,
    npes,
    paths,
    with_breit,
    with_gaunt,
    symmetry="s4",
    tol=1e-8,
    retain_paldus_123=True,
    retain_paldus_456=True,
    only_dg=False,
):
    """Out-of-core binary 2e dump; same record set as ``dump_2e_outcore_opt``.

    ``paths`` maps the Coulomb tags (``LLLL``, ``SSSS``, ``LLSS``, ``SSLL``)
    and the two Breit/Gaunt tags (``LSLS``, ``LSSL``) to the HDF5 files written
    by ``_build_coulomb_outcore`` / ``_build_breit_outcore``.
    """
    import h5py

    _check_index_limit(npes)
    c1 = 0.5 / lib.param.LIGHT_SPEED
    if with_breit:
        sign = 1.0
    elif with_gaunt:
        sign = -1.0
    else:
        sign = 0.0
    write_b = bool(with_breit)
    s4 = symmetry == "s4"
    if symmetry not in ("s4", "s1"):
        raise ValueError("Unknown symmetry %s" % symmetry)

    files = {tag: h5py.File(paths[tag], "r") for tag in paths}
    dsets = {tag: files[tag]["eri_mo"] for tag in files}
    try:
        for i in range(npes):
            row0 = i * npes
            nrows = i + 1 if s4 else npes
            slabs = {tag: np.asarray(dsets[tag][row0 : row0 + nrows]) for tag in dsets}
            for j in range(nrows):
                coul = slabs["LLLL"][j].reshape(npes, npes)
                coul = coul + (c1**2) * slabs["SSLL"][j].reshape(npes, npes)
                coul = coul + (c1**2) * slabs["LLSS"][j].reshape(npes, npes)
                coul = coul + (c1**4) * slabs["SSSS"][j].reshape(npes, npes)
                if write_b:
                    breit = slabs["LSLS"][j].reshape(npes, npes)
                    breit = breit + slabs["LSSL"][j].reshape(npes, npes)
                    ji = j * npes + i
                    breit = breit + _read_row(dsets["LSLS"], ji, npes).T.conj()
                    breit = breit + _read_row(dsets["LSSL"], ji, npes).T.conj()
                    breit = sign * (c1**2) * breit
                else:
                    breit = None
                if s4:
                    coul = coul[: i + 1, :]
                    if breit is not None:
                        breit = breit[: i + 1, :]
                _emit_row_binary(
                    writer,
                    i,
                    j,
                    coul,
                    breit,
                    tol,
                    retain_paldus_123,
                    retain_paldus_456,
                    only_dg,
                    npes,
                )
    finally:
        for h5 in files.values():
            h5.close()
    return writer


# ---------------------------------------------------------------------------
# 1e and ECORE records
# ---------------------------------------------------------------------------


def dump_h1e_binary(writer, h1e, tol=1e-8, only_dg=False):
    """One-electron records (``r = s = 0``), same selection as the ASCII dump."""
    npes = h1e.shape[0]
    _check_index_limit(npes)
    idx = np.arange(npes)
    mask = np.abs(h1e) > tol
    if only_dg:
        mask &= (idx[:, None] // 2) == (idx[None, :] // 2)
    ii, jj = np.nonzero(mask)
    zeros = np.zeros(ii.size, dtype=np.int64)
    writer.add_records(ii + 1, jj + 1, zeros, zeros, h1e[ii, jj])
    return writer


def dump_ecore_binary(writer, nuc):
    """Nuclear repulsion record (``p = q = r = s = 0``)."""
    writer.add_records(
        np.zeros(1, dtype=np.int64),
        np.zeros(1, dtype=np.int64),
        np.zeros(1, dtype=np.int64),
        np.zeros(1, dtype=np.int64),
        np.asarray([nuc], dtype=np.float64),
    )
    return writer


# ---------------------------------------------------------------------------
# main entry points
# ---------------------------------------------------------------------------


def _resolve_orbsym_ID(orbsym_ID, nmo, n2c, with_nes):
    """Same ``ORBSYM`` handling as :func:`FCIDUMP_Rela4C_opt`."""
    if orbsym_ID is None:
        return [0] * nmo
    if len(orbsym_ID) == n2c:
        if with_nes:
            return orbsym_ID[:nmo]
        return orbsym_ID[n2c // 2 : n2c // 2 + nmo]
    if len(orbsym_ID) == (n2c // 2):
        if with_nes:
            raise RuntimeError(
                "orbsym_ID should have length n2c//2 + npes//2 when with_nes=True"
            )
        return orbsym_ID[:nmo]
    raise RuntimeError(
        "orbsym_ID has %d entries, expected n2c=%d or n2c//2=%d"
        % (len(orbsym_ID), n2c, n2c // 2)
    )


def FCIDUMP_Rela4C_binary(
    mol,
    my_RDHF,
    with_breit=None,
    npes=None,
    filename="fcidump",
    binary_filename=None,
    mode="incore",
    orbsym_ID=None,
    IsComplex=True,
    tol=1e-8,
    debug=False,
    no_2e=False,
    retain_paldus_123=True,
    retain_paldus_456=True,
    only_dg=False,
    with_nes=False,
    write_ascii_integrals=False,
    verbose=0,
):
    """Dump four-component integrals: ASCII head + binary integral records.

    Drop-in counterpart of
    :func:`pyscf_util.Relativisitc.integral_4C_opt.FCIDUMP_Rela4C_opt` with the
    same arguments, plus:

    Keyword Args:
        binary_filename (str): binary integral file.  Defaults to
            ``filename + ".bin"``.
        write_ascii_integrals (bool): also write the classic ASCII integral
            records into ``filename`` (default: head only).  Handy to verify
            the binary file, but it doubles the dump time and disk usage.

    Returns:
        ``(int2e_coulomb, int2e_breit)`` when ``debug=True``, else
        ``(None, None)`` -- identical to the ASCII writers.

    Notes:
        * the C++ reader takes the *head* from ``filename`` and *everything
          else* from ``binary_filename``, hence the binary file contains the
          2e, the 1e and the ECORE records;
        * ``IsComplex`` selects the record size: use ``IsComplex=False`` for
          the real ``*_D2H`` integral managers (16-byte records), and
          ``IsComplex=True`` only for managers whose scalar type is complex.
    """
    assert mode in ["original", "incore", "outcore"]

    n2c = mol.nao_2c()
    if npes is None:
        npes = n2c

    mo_spinor, mo_L, mo_S = _mo_parts(mol, my_RDHF, npes, with_nes)
    ndim = mo_spinor.shape[1]  # number of spinors actually dumped
    _check_index_limit(ndim)

    if binary_filename is None:
        binary_filename = filename + ".bin"

    # the head/orbsym must use the same "nmo" convention as the ASCII writers
    nmo = npes // 2
    if with_nes:
        nmo += n2c // 2

    hcore = my_RDHF.get_hcore()
    h1e = mo_spinor.conj().T @ hcore @ mo_spinor

    if with_breit is None:
        with_breit = my_RDHF.with_breit
    with_gaunt = my_RDHF.with_gaunt

    c1 = 0.5 / lib.param.LIGHT_SPEED

    workdir = tempfile.mkdtemp(prefix="rela4c_binary_")
    outcore_paths = None
    int2e_coulomb = None
    int2e_breit = None
    t0 = time.time()
    try:
        if mode == "original":
            # reference/debug path: reuse the original (slow) transformation
            from pyscf_util.Relativisitc.integral_4C import FCIDUMP_Rela4C

            _log(verbose, "  original reference transform ...")
            int2e_coulomb, int2e_breit = FCIDUMP_Rela4C(
                mol,
                my_RDHF,
                with_breit=with_breit,
                npes=npes,
                filename=os.path.join(workdir, "reference_ascii"),
                mode="original",
                orbsym_ID=orbsym_ID,
                IsComplex=IsComplex,
                tol=tol,
                debug=True,
                no_2e=no_2e,
                retain_paldus_123=retain_paldus_123,
                retain_paldus_456=retain_paldus_456,
                only_dg=only_dg,
                with_nes=with_nes,
            )
        elif mode == "incore":
            _log(verbose, "  Coulomb transform ...")
            int2e_coulomb = _build_coulomb_incore(mol, mo_L, mo_S, c1, workdir, verbose)
            if with_breit or with_gaunt:
                terms = _BREIT_TERMS if with_breit else _GAUNT_TERMS
                _log(verbose, "  Breit/Gaunt transform ...")
                int2e_breit = _build_breit_incore(
                    mol, mo_L, mo_S, c1, terms, workdir, verbose
                )
        else:  # outcore
            prefix = os.path.join(workdir, "eri")
            _log(verbose, "  Coulomb transform ...")
            outcore_paths = _build_coulomb_outcore(mol, mo_L, mo_S, prefix, verbose)
            if with_breit or with_gaunt:
                terms = _BREIT_TERMS if with_breit else _GAUNT_TERMS
                _log(verbose, "  Breit/Gaunt transform ...")
                outcore_paths.update(
                    _build_breit_outcore(mol, mo_L, mo_S, terms, prefix, verbose)
                )
        _log(verbose, "  transforms done in %.2fs" % (time.time() - t0))

        orbsym_ID = _resolve_orbsym_ID(orbsym_ID, nmo, n2c, with_nes)
        nelec = mol.nelectron
        ms = 0
        nuc = mol.get_enuc()
        symmetry = "s1" if debug else "s4"

        t_dump = time.time()

        # ---- ASCII file: head (and optionally the integrals) -----------
        with open(filename, "w") as fout:
            tools.fcidump.write_head(fout, nmo, nelec, ms, orbsym_ID)
            if write_ascii_integrals:
                if not no_2e:
                    if mode == "outcore":
                        dump_2e_outcore_opt(
                            fout,
                            ndim,
                            outcore_paths,
                            with_breit,
                            with_gaunt,
                            IsComplex,
                            symmetry=symmetry,
                            tol=tol,
                            retain_paldus_123=retain_paldus_123,
                            retain_paldus_456=retain_paldus_456,
                            only_dg=only_dg,
                        )
                    else:
                        dump_2e_opt(
                            fout,
                            int2e_coulomb,
                            int2e_breit,
                            with_breit,
                            IsComplex,
                            symmetry=symmetry,
                            tol=tol,
                            retain_paldus_123=retain_paldus_123,
                            retain_paldus_456=retain_paldus_456,
                            only_dg=only_dg,
                        )
                _dump_h1e_opt(fout, h1e, tol, IsComplex, only_dg, ndim)
                if IsComplex:
                    fout.write(_FMT_ECORE_C % nuc)
                else:
                    fout.write(_FMT_ECORE_R % nuc)

        # ---- binary file: all integrals --------------------------------
        with BinaryIntegralWriter(binary_filename, is_complex=IsComplex) as writer:
            if not no_2e:
                if mode == "outcore":
                    dump_2e_outcore_binary(
                        writer,
                        ndim,
                        outcore_paths,
                        with_breit,
                        with_gaunt,
                        symmetry=symmetry,
                        tol=tol,
                        retain_paldus_123=retain_paldus_123,
                        retain_paldus_456=retain_paldus_456,
                        only_dg=only_dg,
                    )
                else:
                    dump_2e_binary(
                        writer,
                        int2e_coulomb,
                        int2e_breit,
                        with_breit,
                        symmetry=symmetry,
                        tol=tol,
                        retain_paldus_123=retain_paldus_123,
                        retain_paldus_456=retain_paldus_456,
                        only_dg=only_dg,
                    )
            dump_h1e_binary(writer, h1e, tol=tol, only_dg=only_dg)
            dump_ecore_binary(writer, nuc)
            nrecords = writer.nrecords

        _log(
            verbose,
            "  dump done in %.2fs (%d records, %.1f MB -> %s)"
            % (
                time.time() - t_dump,
                nrecords,
                nrecords * (24 if IsComplex else 16) / 2.0**20,
                binary_filename,
            ),
        )
    finally:
        shutil.rmtree(workdir, ignore_errors=True)

    if debug:
        return int2e_coulomb, int2e_breit
    return None, None


def FCIDUMP_Rela4C_SU2_binary(
    mol,
    my_RDHF,
    with_breit=False,
    npes=None,
    filename="fcidump",
    binary_filename=None,
    mode="incore",
    debug=False,
    no_2e=False,
    with_nes=False,
    tol=1e-8,
    write_ascii_integrals=False,
    verbose=0,
):
    """Jz/parity adapted (D2H) four-component dump in binary format.

    Same orbital adaptation as
    :func:`pyscf_util.Relativisitc.integral_4C_opt.FCIDUMP_Rela4C_SU2_opt`, but
    the integrals go to ``binary_filename`` (default ``filename + ".bin"``)
    and the ASCII file only keeps the head.

    This is the writer matching ``iCIPT2_Breit_D2H_ERICompressed``: the
    real (``IsComplex=False``) records are read by
    ``Integrals4C_D2h_limited_memory`` (``double``), i.e. 16-byte records::

        ./iCIPT2_Breit_D2H_ERICompressed  input  FCIDUMP_F  FCIDUMP_F.bin

    Returns:
        ``(int2e_coulomb, int2e_breit, mo_parity, mo_coeff)``.
    """
    from pyscf_util.Relativisitc.double_group import (
        _atom_Jz_adapted,
        _atom_spinor_spatial_parity,
        atm_d2h_symmetry_adapt_mo_coeff,
    )

    ######## adapt the molecular orbitals to the Jz symmetry and tr symmetry ########

    mo_coeff = _atom_Jz_adapted(mol, my_RDHF.mo_coeff, my_RDHF.mo_energy, debug=debug)
    mo_pes = atm_d2h_symmetry_adapt_mo_coeff(mol, mo_coeff, debug=debug)
    mo_nes = atm_d2h_symmetry_adapt_mo_coeff(mol, mo_coeff, True, debug=debug)
    n2c = mol.nao_2c()
    mo_coeff[:, n2c:] = mo_pes
    mo_coeff[:, :n2c] = mo_nes

    mo_parity = _atom_spinor_spatial_parity(mol, mo_coeff, debug=debug)

    if debug:
        print("mo_parity = ", mo_parity)
    mo_parity = [x for id, x in enumerate(mo_parity) if id % 2 == 0]
    if debug:
        print("mo_parity = ", mo_parity)

    my_RDHF.mo_coeff = mo_coeff

    if npes is None:
        npes = n2c

    if binary_filename is None:
        binary_filename = filename + ".bin"

    coulomb, breit = FCIDUMP_Rela4C_binary(
        mol,
        my_RDHF,
        with_breit=with_breit,
        npes=npes,
        filename=filename,
        binary_filename=binary_filename,
        mode=mode,
        orbsym_ID=mo_parity,
        IsComplex=False,
        tol=tol,
        debug=debug,
        no_2e=no_2e,
        write_ascii_integrals=write_ascii_integrals,
        with_nes=with_nes,
        verbose=verbose,
    )

    if debug:
        # complex-valued cross-check file (24-byte records)
        FCIDUMP_Rela4C_binary(
            mol,
            my_RDHF,
            with_breit=with_breit,
            filename=filename + "_complex",
            binary_filename=filename + "_complex.bin",
            mode=mode,
            orbsym_ID=None,
            IsComplex=True,
            tol=tol,
            debug=debug,
            no_2e=no_2e,
            write_ascii_integrals=write_ascii_integrals,
            npes=npes,
            with_nes=with_nes,
            verbose=verbose,
        )

    return coulomb, breit, mo_parity, mo_coeff


# ---------------------------------------------------------------------------
# readers / verification helpers
# ---------------------------------------------------------------------------


def iter_binary_records(filename, is_complex=False, chunk=1 << 20):
    """Yield the records of a binary integral file as structured arrays."""
    dtype = _record_dtype(is_complex)
    size = os.path.getsize(filename)
    if size % dtype.itemsize != 0:
        raise ValueError(
            "%s: file size %d is not a multiple of the record size %d -- "
            "wrong IsComplex? (real records are 16 bytes, complex 24 bytes)"
            % (filename, size, dtype.itemsize)
        )
    total = size // dtype.itemsize
    with open(filename, "rb") as fin:
        remaining = total
        while remaining > 0:
            count = int(min(chunk, remaining))
            data = np.fromfile(fin, dtype=dtype, count=count)
            if data.size != count:
                raise ValueError("%s: truncated binary integral file" % filename)
            remaining -= count
            yield data


def read_binary_records(filename, is_complex=False):
    """Read a whole binary integral file (structured numpy array)."""
    dtype = _record_dtype(is_complex)
    size = os.path.getsize(filename)
    if size % dtype.itemsize != 0:
        raise ValueError(
            "%s: file size %d is not a multiple of the record size %d -- "
            "wrong IsComplex? (real records are 16 bytes, complex 24 bytes)"
            % (filename, size, dtype.itemsize)
        )
    return np.fromfile(filename, dtype=dtype)


def binary_summary(filename, is_complex=False):
    """Count the record types of a binary integral file."""
    dtype = _record_dtype(is_complex)
    size = os.path.getsize(filename)
    if size % dtype.itemsize != 0:
        raise ValueError(
            "%s: file size %d is not a multiple of the record size %d"
            % (filename, size, dtype.itemsize)
        )
    counts = {
        "nrecords": size // dtype.itemsize,
        "record_size": dtype.itemsize,
        "n_2e": 0,
        "n_1e": 0,
        "n_ecore": 0,
        "bytes": size,
    }
    for records in iter_binary_records(filename, is_complex):
        r = records["r"]
        s_ = records["s"]
        p = records["p"]
        q = records["q"]
        two_e = r != 0
        one_e = (r == 0) & (s_ == 0) & (p != 0)
        ecore = (p == 0) & (q == 0) & (r == 0) & (s_ == 0)
        counts["n_2e"] += int(np.count_nonzero(two_e))
        counts["n_1e"] += int(np.count_nonzero(one_e))
        counts["n_ecore"] += int(np.count_nonzero(ecore))
    return counts


def binary_to_fcidump(
    binary_filename,
    outname,
    is_complex=False,
    head_filename=None,
    norb=None,
    nelec=None,
    ms2=0,
    orbsym=None,
):
    """Convert a binary integral file back to a classic ASCII FCIDUMP.

    The head is either copied from ``head_filename`` (an existing FCIDUMP, e.g.
    the head-only file written by this module) or generated from
    ``norb`` / ``nelec`` / ``orbsym``.
    """
    dtype = _record_dtype(is_complex)
    size = os.path.getsize(binary_filename)
    if size % dtype.itemsize != 0:
        raise ValueError(
            "%s: file size %d is not a multiple of the record size %d -- "
            "wrong IsComplex? (real records are 16 bytes, complex 24 bytes)"
            % (binary_filename, size, dtype.itemsize)
        )
    nrecords = size // dtype.itemsize

    with open(outname, "w") as fout:
        if head_filename is not None:
            with open(head_filename) as fin:
                for line in fin:
                    fout.write(line)
                    if "&END" in line:
                        break
        else:
            if norb is None or nelec is None:
                raise ValueError(
                    "binary_to_fcidump needs either head_filename or "
                    "norb/nelec/orbsym"
                )
            tools.fcidump.write_head(fout, int(norb), int(nelec), int(ms2), orbsym)

        chunk = 1 << 16
        for records in iter_binary_records(binary_filename, is_complex, chunk=chunk):
            p = records["p"].tolist()
            q = records["q"].tolist()
            r = records["r"].tolist()
            s = records["s"].tolist()
            if is_complex:
                re = records["vr"].tolist()
                im = records["vi"].tolist()
                fmt = _ASCII_FMT_COMPLEX
                lines = [
                    fmt % (re[i], im[i], p[i], q[i], r[i], s[i]) for i in range(len(p))
                ]
            else:
                values = records["v"].tolist()
                fmt = _ASCII_FMT_REAL
                lines = [
                    fmt % (values[i], p[i], q[i], r[i], s[i]) for i in range(len(p))
                ]
            fout.write("".join(lines))
    return {"filename": outname, "nrecords": int(nrecords)}


def fcidump_to_binary(
    ascii_filename,
    head_filename,
    binary_filename,
    is_complex=None,
    verbose=0,
):
    """Split a classic ASCII FCIDUMP into (head-only ASCII, binary records).

    This is the inverse of :func:`binary_to_fcidump`: it lets an *existing*
    four-component FCIDUMP be reused by the binary readers without recomputing
    the integrals.  The record order and the values are preserved.

    Args:
        ascii_filename: classic FCIDUMP (head + ASCII integral records).
        head_filename: output file that receives the head only.
        binary_filename: output binary file.

    Kwargs:
        is_complex: ``None`` (default) auto-detects from the first record
            (6 columns = complex, 5 columns = real).

    Returns:
        ``{"head": ..., "binary": ..., "nrecords": n, "is_complex": bool}``
    """
    nrecords = 0
    head_done = False
    first_record = None

    with open(ascii_filename) as fin, open(head_filename, "w") as fhead:
        for line in fin:
            fhead.write(line)
            if "&END" in line:
                head_done = True
                break
        if not head_done:
            raise ValueError("%s: no '&END' in the head" % ascii_filename)

        for line in fin:
            if line.strip():
                first_record = line
                break
        if first_record is not None:
            n_columns = len(first_record.split())
            detected = n_columns == 6
            if is_complex is None:
                is_complex = detected
            elif detected != bool(is_complex):
                raise ValueError(
                    "%s: first record has %d columns but is_complex=%s"
                    % (ascii_filename, n_columns, is_complex)
                )
        elif is_complex is None:
            is_complex = False

        with BinaryIntegralWriter(
            binary_filename, is_complex=bool(is_complex)
        ) as writer:
            chunk = []

            def _flush():
                nonlocal chunk, nrecords
                if not chunk:
                    return
                values = np.asarray([row[0] for row in chunk])
                indices = np.asarray([row[1:] for row in chunk], dtype=np.int64)
                writer.add_records(
                    indices[:, 0],
                    indices[:, 1],
                    indices[:, 2],
                    indices[:, 3],
                    values,
                )
                nrecords += len(chunk)
                chunk = []

            def _consume(line):
                nonlocal chunk
                tokens = line.split()
                if not tokens:
                    return
                if is_complex:
                    if len(tokens) != 6:
                        raise ValueError(
                            "%s: expected 6 columns, got %r" % (ascii_filename, line)
                        )
                    value = complex(float(tokens[0]), float(tokens[1]))
                    p, q, r, s = (int(x) for x in tokens[2:6])
                else:
                    if len(tokens) != 5:
                        raise ValueError(
                            "%s: expected 5 columns, got %r" % (ascii_filename, line)
                        )
                    value = float(tokens[0])
                    p, q, r, s = (int(x) for x in tokens[1:5])
                chunk.append((value, p, q, r, s))
                if len(chunk) >= (1 << 16):
                    _flush()

            if first_record is not None:
                _consume(first_record)
            for line in fin:
                _consume(line)
            _flush()

    _log(
        verbose,
        "converted %d records (%s) -> %s + %s"
        % (
            nrecords,
            "complex" if is_complex else "real",
            head_filename,
            binary_filename,
        ),
    )
    return {
        "head": head_filename,
        "binary": binary_filename,
        "nrecords": int(nrecords),
        "is_complex": bool(is_complex),
    }


def parse_fcidump_records(path):
    """Parse an ASCII FCIDUMP into ``(records, h1e, ecore, head)``.

    ``records`` maps ``(p, q, r, s)`` (1-based, 2e only) to a complex value;
    ``h1e`` maps ``(p, q)`` to a complex value; ``head`` holds the raw head
    text up to and including ``&END``.
    """
    import re

    with open(path) as fin:
        lines = fin.read().splitlines()

    end = None
    for idx, line in enumerate(lines):
        if "&END" in line:
            end = idx
            break
    if end is None:
        raise ValueError("%s: no '&END' in the head" % path)

    records = {}
    h1e = {}
    ecore = None
    for line in lines[end + 1 :]:
        tokens = line.split()
        if len(tokens) < 5:
            continue
        if len(tokens) == 6:
            value = complex(float(tokens[0]), float(tokens[1]))
            p, q, r, s = (int(x) for x in tokens[2:6])
        else:
            value = complex(float(tokens[0]), 0.0)
            p, q, r, s = (int(x) for x in tokens[1:5])
        if r == 0 and s == 0:
            if p == 0 and q == 0:
                ecore = value
            else:
                h1e[(p, q)] = value
            continue
        records[(p, q, r, s)] = value
    return records, h1e, ecore, "\n".join(lines[: end + 1])


def compare_binary_with_fcidump(
    binary_filename,
    ascii_filename,
    is_complex=False,
    tol=1e-10,
    verbose=False,
):
    """Compare a binary integral file with the equivalent ASCII FCIDUMP.

    Returns a dictionary with the record counts and the maximum absolute
    deviation of the 2e / 1e / ECORE values; raises ``AssertionError`` if the
    record sets differ or a value is off by more than ``tol``.
    """
    records, h1e_ref, ecore_ref, _head = parse_fcidump_records(ascii_filename)
    data = read_binary_records(binary_filename, is_complex)

    if is_complex:
        values = data["vr"].astype(np.complex128) + 1j * data["vi"]
    else:
        values = data["v"].astype(np.complex128)
    keys = [
        (int(p), int(q), int(r), int(s))
        for p, q, r, s in zip(data["p"], data["q"], data["r"], data["s"])
    ]

    got_2e = {}
    got_1e = {}
    got_ecore = None
    for key, value in zip(keys, values.tolist()):
        p, q, r, s = key
        if r == 0 and s == 0:
            if p == 0 and q == 0:
                got_ecore = value
            else:
                got_1e[(p, q)] = value
            continue
        got_2e[key] = value

    if set(got_2e) != set(records):
        missing = sorted(set(records) - set(got_2e))[:5]
        extra = sorted(set(got_2e) - set(records))[:5]
        raise AssertionError(
            "2e record sets differ: %d in ASCII, %d in binary "
            "(e.g. missing %s, extra %s)" % (len(records), len(got_2e), missing, extra)
        )
    if set(got_1e) != set(h1e_ref):
        raise AssertionError(
            "1e record sets differ: %d in ASCII, %d in binary"
            % (len(h1e_ref), len(got_1e))
        )

    delta_2e = max((abs(got_2e[k] - records[k]) for k in records), default=0.0)
    delta_1e = max((abs(got_1e[k] - h1e_ref[k]) for k in h1e_ref), default=0.0)
    if ecore_ref is None and got_ecore is None:
        delta_ecore = 0.0
    elif ecore_ref is None or got_ecore is None:
        raise AssertionError("ECORE record missing in ASCII or in binary file")
    else:
        delta_ecore = abs(got_ecore - ecore_ref)

    result = {
        "n_2e": len(records),
        "n_1e": len(h1e_ref),
        "delta_2e": delta_2e,
        "delta_1e": delta_1e,
        "delta_ecore": delta_ecore,
    }
    if verbose:
        print(
            "  %-28s %s"
            % ("records (2e/1e)", "%d / %d" % (result["n_2e"], result["n_1e"]))
        )
        print("  %-28s %.3e" % ("max |d 2e|", delta_2e))
        print("  %-28s %.3e" % ("max |d 1e|", delta_1e))
        print("  %-28s %.3e" % ("|d ECORE|", delta_ecore))

    assert delta_2e <= tol, "2e deviation %.3e > %.3e" % (delta_2e, tol)
    assert delta_1e <= tol, "1e deviation %.3e > %.3e" % (delta_1e, tol)
    assert delta_ecore <= tol, "ECORE deviation %.3e > %.3e" % (delta_ecore, tol)
    return result


# ---------------------------------------------------------------------------
# small demo
# ---------------------------------------------------------------------------


def _demo():
    import tempfile as _tempfile

    from pyscf import gto, scf

    mol = gto.M(
        atom="O 0 0 0; H 0 -0.757 0.587; H 0 0.757 0.587",
        basis="3-21g",
        verbose=0,
        charge=0,
        spin=0,
        symmetry="c2v",
    )
    mol.build()
    mf = scf.dhf.RDHF(mol)
    mf.conv_tol = 1e-10
    mf.kernel()
    mf.with_breit = True

    workdir = _tempfile.mkdtemp(prefix="rela4c_binary_demo_")
    head_file = os.path.join(workdir, "FCIDUMP_demo")
    ascii_file = os.path.join(workdir, "FCIDUMP_demo_ascii")
    binary_file = head_file + ".bin"

    # reference ASCII dump (classic records) and the binary dump; the D2H
    # iCIPT2 integral manager is real -> IsComplex=False (16-byte records)
    from pyscf_util.Relativisitc.integral_4C_opt import FCIDUMP_Rela4C_opt

    FCIDUMP_Rela4C_opt(
        mol, mf, with_breit=True, npes=20, filename=ascii_file, IsComplex=False
    )
    FCIDUMP_Rela4C_binary(
        mol,
        mf,
        with_breit=True,
        npes=20,
        filename=head_file,
        mode="incore",
        IsComplex=False,
        write_ascii_integrals=False,
        verbose=1,
    )

    print("binary summary:", binary_summary(binary_file, is_complex=False))
    print(compare_binary_with_fcidump(binary_file, ascii_file, verbose=True))
    print(
        binary_to_fcidump(
            binary_file, head_file + ".roundtrip", head_filename=head_file
        )
    )

    with open(head_file) as fin:
        print("".join(fin.readlines()))
    return binary_file


if __name__ == "__main__":
    _demo()
