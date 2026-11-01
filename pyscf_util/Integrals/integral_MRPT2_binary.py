#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Binary MRPT2 integral dump (FCIDUMP-style ASCII head + binary integrals).

Motivation
----------
FCIDUMP files written by the MRPT2 writers are dominated by the two-electron
integrals (a single `` %.16g %4d %4d %4d %4d`` record occupies 30+ bytes).
The binary format implemented here keeps the *head* in ASCII -- exactly like
FCIDUMP, so that ``norb``, ``nelec``, ``ms2``, ``orbsym`` etc. remain human
readable and tool friendly -- and stores the integrals as raw ``float64``
values with ``uint16`` orbital indices (16 bytes per record).

File layout
-----------
::

    [ASCII head]                      # comma separated key=value tokens
     &FCI NORB=%4d,NELEC=%2d,MS2=%d,
      ORBSYM=<nmo integers>
      ISYM=1,
      NFZC=%4d,NACT=%4d,NVIR=%4d,TOL=%.17g,
      ECORE=%.17g,FMT=MRPT2BIN1,IDXBASE=0,MODE=%s,
     &END
    [32 byte preamble]                # little endian
     8s  b"MRPT2BIN"
     u64 format version (1)
     u64 nmo
     u64 reserved (0)
    [h1e]        (nmo*(nmo+1)/2) x "<f8", lower triangle, row major (i >= j)
    [ERI]        nint x 16 bytes:
                 u16 p | u16 q | u16 r | u16 s | f8 value
    [16 byte trailer]
     8s  b"MRPT2END"
     u64 nint

Conventions
-----------
* indices are 0-based (``IDXBASE=0`` in the head);
* each record is in canonical 8-fold form: ``p >= q``, ``r >= s`` and
  ``combine2(p, q) >= combine2(r, s)`` (the same packing as
  :func:`pyscf_util.misc.misc._combine4`), with ``value = (p q | r s)``;
* 8-fold keys are unique inside a file;
* the order of the records is *not* part of the format contract (the sparse
  writers emit ascending keys, the fast writer emits block order).  Readers
  must not rely on the order;
* ``nmo`` must be ``<= 65536`` (``uint16`` index range);
* little endian only.

Public API
----------
Writers::

    write_mrpt2_binary(filename, h1e, chunks, *, nmo, nelec, ...)
    dump_mrpt2_binary(mol, scf, mo_coeff, nfzc, nact, nvir, ...)          # incore
    dump_mrpt2_binary_outcore(mol, scf, mo_coeff, nfzc, nact, nvir, ...)  # outcore
    dump_mrpt2_binary_fast(mol, scf, mo_coeff, nfzc, nact, nvir, ...)     # block set

Readers::

    read_mrpt2_binary_header(filename)
    iter_mrpt2_binary(filename, chunk=1 << 20)
    load_mrpt2_binary(filename)
    mrpt2_binary_to_fcidump(filename, outname)
    mrpt2_binary_to_pyscf_dict(filename)   # same structure as fcidump.read()

Record set parity with the ASCII writers
----------------------------------------
* ``dump_mrpt2_binary(..., mode="mrpt2")``     == ``fcidump_mrpt2``
* ``dump_mrpt2_binary(..., mode="full")``      == all non-zero 8-fold integrals
* ``dump_mrpt2_binary_outcore(...)``           == ``fcidump_mrpt2_outcore``
* ``dump_mrpt2_binary_fast(...)``              == ``fcidump_mrpt2_incore_fast``
"""

from __future__ import annotations

import os
import re
import struct
import tempfile

import h5py
import numpy as np
from pyscf import ao2mo, tools

from pyscf_util.Integrals.integral_MRPT2_opt import (
    _outcore_lookup,
    _prepare,
    fast_chunks,
    full_chunks,
    incore_lookup,
    mask_chunks,
    orbsym_ids,
    prepare_erifile,
)

__all__ = [
    "MRPT2BIN_MAGIC",
    "MRPT2BIN_TRAILER",
    "MRPT2BIN_VERSION",
    "MRPT2_ERI_DTYPE",
    "write_mrpt2_binary",
    "dump_mrpt2_binary",
    "dump_mrpt2_binary_outcore",
    "dump_mrpt2_binary_fast",
    "read_mrpt2_binary_header",
    "iter_mrpt2_binary",
    "load_mrpt2_binary",
    "mrpt2_binary_to_fcidump",
    "mrpt2_binary_to_pyscf_dict",
]

MRPT2BIN_MAGIC = b"MRPT2BIN"
MRPT2BIN_TRAILER = b"MRPT2END"
MRPT2BIN_VERSION = 1
MRPT2BIN_FMT = "MRPT2BIN1"
MRPT2_MAX_NMO = 1 << 16  # uint16: indices 0 .. 65535

MRPT2_ERI_DTYPE = np.dtype(
    [("p", "<u2"), ("q", "<u2"), ("r", "<u2"), ("s", "<u2"), ("v", "<f8")]
)
assert MRPT2_ERI_DTYPE.itemsize == 16

_FLOAT_FORMAT = tools.fcidump.DEFAULT_FLOAT_FORMAT
_HEAD_PREAMBLE = struct.Struct("<8sQQQ")  # 32 bytes
_HEAD_TRAILER = struct.Struct("<8sQ")  # 16 bytes
_RECORDS_PER_FLUSH = 1 << 20


# ---------------------------------------------------------------------------
# ASCII head
# ---------------------------------------------------------------------------


def _head_text(nmo, nelec, ms2, orbsym, nfzc, nact, nvir, tol, ecore, mode):
    if orbsym is None:
        orbsym = [1] * nmo
    if len(orbsym) != nmo:
        raise ValueError(
            "orbsym has %d entries, expected nmo=%d" % (len(orbsym), nmo)
        )
    return (
        " &FCI NORB=%4d,NELEC=%2d,MS2=%d,\n" % (nmo, nelec, ms2)
        + "  ORBSYM=%s\n" % ",".join(str(int(x)) for x in orbsym)
        + "  ISYM=1,\n"
        + "  NFZC=%4d,NACT=%4d,NVIR=%4d,TOL=%.17g,\n"
        % (nfzc, nact, nvir, tol)
        + "  ECORE=%.17g,FMT=%s,IDXBASE=0,MODE=%s,\n" % (ecore, MRPT2BIN_FMT, mode)
        + " &END\n"
    )


def _find_int(text, key, default=None):
    match = re.search(r"\b%s\s*=\s*(-?\d+)" % re.escape(key), text)
    if match is None:
        return default
    return int(match.group(1))


def _find_float(text, key, default=None):
    pattern = r"\b%s\s*=\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)" % re.escape(key)
    match = re.search(pattern, text)
    if match is None:
        return default
    return float(match.group(1))


def _find_str(text, key, default=None):
    match = re.search(r"\b%s\s*=\s*([A-Za-z0-9_.+-]+)" % re.escape(key), text)
    if match is None:
        return default
    return match.group(1)


def _parse_head(text):
    meta = {}
    meta["nmo"] = _find_int(text, "NORB")
    meta["nelec"] = _find_int(text, "NELEC")
    meta["ms2"] = _find_int(text, "MS2", 0)
    meta["isym"] = _find_int(text, "ISYM", 1)
    meta["nfzc"] = _find_int(text, "NFZC", 0)
    meta["nact"] = _find_int(text, "NACT", 0)
    meta["nvir"] = _find_int(text, "NVIR", 0)
    meta["tol"] = _find_float(text, "TOL", 0.0)
    meta["ecore"] = _find_float(text, "ECORE", 0.0)
    meta["fmt"] = _find_str(text, "FMT", "")
    meta["idxbase"] = _find_int(text, "IDXBASE", 0)
    meta["mode"] = _find_str(text, "MODE", "")
    match = re.search(r"\bORBSYM\s*=\s*([0-9, \t]+)", text)
    if match is None:
        meta["orbsym"] = None
    else:
        meta["orbsym"] = [
            int(x) for x in re.split(r"[,\s]+", match.group(1).strip()) if x
        ]
    return meta


# ---------------------------------------------------------------------------
# writer
# ---------------------------------------------------------------------------


def write_mrpt2_binary(
    filename,
    h1e,
    chunks,
    *,
    nmo,
    nelec,
    ms2=0,
    orbsym=None,
    nfzc=0,
    nact=0,
    nvir=0,
    tol=0.0,
    ecore=0.0,
    mode="mrpt2",
    verbose=0,
):
    """Write one binary MRPT2 integral file.

    Parameters
    ----------
    filename : str
    h1e : (nmo, nmo) ndarray
        One-electron Hamiltonian in the MO basis (only the lower triangle is
        stored).
    chunks : iterable of ``(i, j, k, l, vals)``
        0-based index arrays and integral values; see
        :func:`pyscf_util.Integrals.integral_MRPT2_opt.mask_chunks`.
    """
    if not isinstance(mode, str) or not mode:
        raise ValueError("mode must be a non-empty string")
    if nmo < 0:
        raise ValueError("nmo must be non-negative")
    if nmo > MRPT2_MAX_NMO:
        raise ValueError(
            "nmo=%d cannot be stored with uint16 orbital indices (max %d)"
            % (nmo, MRPT2_MAX_NMO)
        )
    h1e = np.ascontiguousarray(h1e, dtype=np.float64)
    if h1e.shape != (nmo, nmo):
        raise ValueError("h1e must have shape (%d, %d), got %s" % (nmo, nmo, h1e.shape))

    head = _head_text(nmo, nelec, ms2, orbsym, nfzc, nact, nvir, tol, ecore, mode)
    head_bytes = head.encode("ascii")
    lower = np.ascontiguousarray(h1e[np.tril_indices(nmo)], dtype="<f8")

    nint = 0
    with open(filename, "wb") as fout:
        fout.write(head_bytes)
        fout.write(
            _HEAD_PREAMBLE.pack(
                MRPT2BIN_MAGIC, MRPT2BIN_VERSION, int(nmo), 0
            )
        )
        fout.write(lower.tobytes())
        eri_offset = fout.tell()

        buffer = []
        buffered = 0

        def flush():
            nonlocal buffer, buffered
            if not buffer:
                return
            fout.write(np.concatenate(buffer).tobytes())
            buffer = []
            buffered = 0

        for i, j, k, l, vals in chunks:
            vals = np.asarray(vals, dtype=np.float64)
            count = int(vals.size)
            if count == 0:
                continue
            indices = (i, j, k, l)
            for arr in indices:
                arr = np.asarray(arr)
                if arr.size != count:
                    raise ValueError(
                        "index and value arrays must have the same length"
                    )
                if arr.size and (arr.min() < 0 or arr.max() >= nmo):
                    raise ValueError(
                        "orbital index out of range [0, %d): %d"
                        % (nmo, int(arr.max()))
                    )
            record = np.empty(count, dtype=MRPT2_ERI_DTYPE)
            record["p"] = np.asarray(i, dtype=np.uint16)
            record["q"] = np.asarray(j, dtype=np.uint16)
            record["r"] = np.asarray(k, dtype=np.uint16)
            record["s"] = np.asarray(l, dtype=np.uint16)
            record["v"] = vals
            buffer.append(record)
            buffered += count
            nint += count
            if buffered >= _RECORDS_PER_FLUSH:
                flush()
        flush()
        fout.write(_HEAD_TRAILER.pack(MRPT2BIN_TRAILER, int(nint)))

    if verbose > 0:
        print(
            "wrote %d binary records to %s (%.1f MB)"
            % (nint, filename, os.path.getsize(filename) / 2.0**20)
        )
    return {
        "filename": filename,
        "nint": int(nint),
        "bytes": int(os.path.getsize(filename)),
        "eri_offset": int(eri_offset),
    }


# ---------------------------------------------------------------------------
# reader
# ---------------------------------------------------------------------------


def read_mrpt2_binary_header(filename):
    """Parse and validate the head + trailer; return offsets and metadata."""
    size = os.path.getsize(filename)
    with open(filename, "rb") as fin:
        head_lines = []
        while True:
            line = fin.readline()
            if not line:
                raise ValueError(
                    "%s: invalid MRPT2 binary file (no '&END' in the head)" % filename
                )
            try:
                text = line.decode("ascii")
            except UnicodeDecodeError:
                raise ValueError(
                    "%s: invalid MRPT2 binary head (non ASCII content)" % filename
                )
            head_lines.append(text)
            if "&END" in text.upper():
                break
            if len(head_lines) > 64:
                raise ValueError(
                    "%s: invalid MRPT2 binary head ('&END' not found)" % filename
                )
        head_end = fin.tell()

        meta = _parse_head("".join(head_lines))
        if meta["fmt"] != MRPT2BIN_FMT:
            raise ValueError(
                "%s: unsupported format %r (expected %r)"
                % (filename, meta["fmt"], MRPT2BIN_FMT)
            )
        if meta["idxbase"] != 0:
            raise ValueError(
                "%s: unsupported IDXBASE=%d (only 0-based files are readable)"
                % (filename, meta["idxbase"])
            )
        preamble = fin.read(_HEAD_PREAMBLE.size)
        if len(preamble) != _HEAD_PREAMBLE.size:
            raise ValueError("%s: truncated binary preamble" % filename)
        magic, version, nmo, reserved = _HEAD_PREAMBLE.unpack(preamble)
        if magic != MRPT2BIN_MAGIC:
            raise ValueError(
                "%s: bad magic %r (expected %r)" % (filename, magic, MRPT2BIN_MAGIC)
            )
        if version != MRPT2BIN_VERSION:
            raise ValueError(
                "%s: unsupported format version %d (expected %d)"
                % (filename, version, MRPT2BIN_VERSION)
            )
        if meta["nmo"] is not None and int(meta["nmo"]) != int(nmo):
            raise ValueError(
                "%s: NORB=%d in the head disagrees with nmo=%d in the preamble"
                % (filename, meta["nmo"], nmo)
            )
        nh1 = int(nmo) * (int(nmo) + 1) // 2
        h1_offset = head_end + _HEAD_PREAMBLE.size
        eri_offset = h1_offset + 8 * nh1
        if size < eri_offset + _HEAD_TRAILER.size:
            raise ValueError("%s: truncated file" % filename)

        fin.seek(size - _HEAD_TRAILER.size)
        trailer = fin.read(_HEAD_TRAILER.size)
        tmagic, nint = _HEAD_TRAILER.unpack(trailer)
        if tmagic != MRPT2BIN_TRAILER:
            raise ValueError(
                "%s: bad trailer %r (expected %r)"
                % (filename, tmagic, MRPT2BIN_TRAILER)
            )
        expected = size - _HEAD_TRAILER.size - eri_offset
        if expected != 16 * nint:
            raise ValueError(
                "%s: record count mismatch (trailer says %d records, "
                "payload holds %d bytes)" % (filename, nint, expected)
            )

    meta.update(
        {
            "format_version": int(version),
            "nmo": int(nmo),
            "nint": int(nint),
            "head_end": int(head_end),
            "h1_offset": int(h1_offset),
            "h1_count": int(nh1),
            "eri_offset": int(eri_offset),
            "filesize": int(size),
        }
    )
    return meta


def iter_mrpt2_binary(filename, chunk=1 << 20):
    """Yield the ERI records as structured arrays (see ``MRPT2_ERI_DTYPE``)."""
    header = read_mrpt2_binary_header(filename)
    with open(filename, "rb") as fin:
        fin.seek(header["eri_offset"])
        remaining = header["nint"]
        while remaining > 0:
            count = min(int(chunk), remaining)
            data = np.fromfile(fin, dtype=MRPT2_ERI_DTYPE, count=count)
            if data.size != count:
                raise ValueError("%s: truncated ERI section" % filename)
            remaining -= count
            yield data


def load_mrpt2_binary(filename):
    """Return the full content of a binary MRPT2 file as a dictionary.

    Keys: ``nmo, nelec, ms2, orbsym, isym, nfzc, nact, nvir, tol, ecore,
    mode, h1e, nuc, eri, nint``.  ``eri`` is a structured array with fields
    ``p, q, r, s, v`` (0-based, canonical 8-fold form).
    """
    header = read_mrpt2_binary_header(filename)
    nmo = header["nmo"]
    with open(filename, "rb") as fin:
        fin.seek(header["h1_offset"])
        triangle = np.fromfile(fin, dtype="<f8", count=header["h1_count"])
    if triangle.size != header["h1_count"]:
        raise ValueError("%s: truncated h1e section" % filename)
    h1e = np.zeros((nmo, nmo))
    rows, cols = np.tril_indices(nmo)
    h1e[rows, cols] = triangle
    h1e[cols, rows] = triangle

    if header["nint"]:
        records = np.concatenate(list(iter_mrpt2_binary(filename)))
    else:
        records = np.empty(0, dtype=MRPT2_ERI_DTYPE)

    header = dict(header)
    header.update({"h1e": h1e, "nuc": header["ecore"], "eri": records})
    return header


def mrpt2_binary_to_fcidump(filename, outname):
    """Convert a binary MRPT2 file back to a standard ASCII FCIDUMP."""
    data = load_mrpt2_binary(filename)
    nmo = data["nmo"]
    records = data["eri"]
    with open(outname, "w") as fout:
        tools.fcidump.write_head(
            fout, nmo, data["nelec"], data["ms2"], data["orbsym"]
        )
        fmt = _FLOAT_FORMAT + " %4d %4d %4d %4d\n"
        lines = []
        for start in range(0, len(records), 1 << 15):
            stop = min(len(records), start + (1 << 15))
            chunk = records[start:stop]
            lines = [
                fmt
                % (
                    v,
                    p + 1,
                    q + 1,
                    r + 1,
                    s + 1,
                )
                for v, p, q, r, s in zip(
                    chunk["v"].tolist(),
                    chunk["p"].tolist(),
                    chunk["q"].tolist(),
                    chunk["r"].tolist(),
                    chunk["s"].tolist(),
                )
            ]
            fout.write("".join(lines))
        tools.fcidump.write_hcore(
            fout, data["h1e"], nmo, tol=0.0, float_format=_FLOAT_FORMAT
        )
        fout.write((_FLOAT_FORMAT + "  0  0  0  0\n") % data["ecore"])
    return {"filename": outname, "nrecords": int(len(records))}


def mrpt2_binary_to_pyscf_dict(filename):
    """Return the same dictionary that :func:`pyscf.tools.fcidump.read` returns.

    ``H2`` uses the 8-fold packed layout expected by PySCF (see
    ``pyscf.tools.fcidump.read``), so the result can be handed to any code
    working with ``fcidump.read`` output.
    """
    data = load_mrpt2_binary(filename)
    nmo = data["nmo"]
    records = data["eri"]
    npair = nmo * (nmo + 1) // 2
    h2 = np.zeros(npair * (npair + 1) // 2)
    if len(records):
        ij = records["p"].astype(np.int64) * (records["p"].astype(np.int64) + 1) // 2
        ij += records["q"]
        kl = records["r"].astype(np.int64) * (records["r"].astype(np.int64) + 1) // 2
        kl += records["s"]
        lo = np.minimum(ij, kl)
        hi = np.maximum(ij, kl)
        h2[hi * (hi + 1) // 2 + lo] = records["v"]
    return {
        "H1": data["h1e"],
        "H2": h2,
        "ECORE": data["ecore"],
        "NORB": nmo,
        "NELEC": data["nelec"],
        "MS2": data["ms2"],
        "ORBSYM": data["orbsym"],
        "ISYM": data["isym"],
    }


# ---------------------------------------------------------------------------
# high level dumps
# ---------------------------------------------------------------------------


def _header_kwargs(mol, scf, mo_coeff, nfzc, nact, nvir, tol, mode):
    nmo, mo, h1e, nelec, nuc = _prepare(mol, scf, mo_coeff, nfzc, nact, nvir)
    return (
        dict(
            nmo=nmo,
            nelec=nelec,
            ms2=0,
            orbsym=orbsym_ids(mol, mo),
            nfzc=nfzc,
            nact=nact,
            nvir=nvir,
            tol=tol,
            ecore=nuc,
            mode=mode,
        ),
        mo,
        h1e,
    )


def dump_mrpt2_binary(
    mol,
    scf,
    mo_coeff,
    nfzc,
    nact,
    nvir,
    filename="FCIDUMP.bin",
    tol=1e-10,
    mode="mrpt2",
    verbose=0,
):
    """In-core binary dump; ``mode="mrpt2"`` matches ``fcidump_mrpt2``."""
    if mode not in ("mrpt2", "full"):
        raise ValueError("mode must be 'mrpt2' or 'full', got %r" % (mode,))
    kwargs, mo, h1e = _header_kwargs(mol, scf, mo_coeff, nfzc, nact, nvir, tol, mode)
    nmo = kwargs["nmo"]
    if nmo > MRPT2_MAX_NMO:
        raise ValueError("nmo=%d exceeds the uint16 limit" % nmo)
    int2e4 = ao2mo.full(eri_or_mol=mol, mo_coeff=mo, aosym="4")
    int2e8 = ao2mo.restore(8, int2e4, nmo)
    del int2e4
    lookup = incore_lookup(int2e8)
    if mode == "full":
        chunks = full_chunks(nmo, tol, lookup)
    else:
        chunks = mask_chunks(nmo, nfzc, nact, tol, lookup)
    stats = write_mrpt2_binary(filename, h1e, chunks, verbose=verbose, **kwargs)
    return stats


def dump_mrpt2_binary_outcore(
    mol,
    scf,
    mo_coeff,
    nfzc,
    nact,
    nvir,
    filename="FCIDUMP.bin",
    tol=1e-8,
    mode="mrpt2",
    erifile=None,
    keep_erifile=False,
    verbose=0,
):
    """Out of core binary dump; ``mode="mrpt2"`` matches ``fcidump_mrpt2_outcore``.

    ``erifile`` may reference an existing HDF5 file with the 4-fold MO
    integrals in dataset ``eri_mo`` (reused instead of recomputed).
    """
    if mode not in ("mrpt2", "full"):
        raise ValueError("mode must be 'mrpt2' or 'full', got %r" % (mode,))
    kwargs, mo, h1e = _header_kwargs(mol, scf, mo_coeff, nfzc, nact, nvir, tol, mode)
    nmo = kwargs["nmo"]
    if nmo > MRPT2_MAX_NMO:
        raise ValueError("nmo=%d exceeds the uint16 limit" % nmo)

    erifile, created = prepare_erifile(mol, mo, nmo, erifile, verbose=verbose)
    try:
        with h5py.File(erifile, "r") as h5:
            dset = h5["eri_mo"]
            lookup = _outcore_lookup(dset, nmo)
            if mode == "full":
                chunks = full_chunks(nmo, tol, lookup)
            else:
                chunks = mask_chunks(nmo, nfzc, nact, tol, lookup)
            stats = write_mrpt2_binary(
                filename, h1e, chunks, verbose=verbose, **kwargs
            )
    finally:
        if created and not keep_erifile:
            try:
                os.unlink(erifile)
            except OSError:
                pass
    return stats


def dump_mrpt2_binary_fast(
    mol,
    scf,
    mo_coeff,
    nfzc,
    nact,
    nvir,
    filename="FCIDUMP.bin",
    tol=1e-10,
    full_integrals=False,
    dedupe=True,
    aosym="s4",
    verbose=0,
):
    """Binary dump of the ``fcidump_mrpt2_incore_fast`` integral set."""
    kwargs, mo, h1e = _header_kwargs(mol, scf, mo_coeff, nfzc, nact, nvir, tol, "fast")
    nmo = kwargs["nmo"]
    if nmo > MRPT2_MAX_NMO:
        raise ValueError("nmo=%d exceeds the uint16 limit" % nmo)
    chunks = fast_chunks(
        mol,
        mo,
        nfzc,
        nact,
        nvir,
        tol,
        full_integrals=full_integrals,
        dedupe=dedupe,
        aosym=aosym,
        verbose=verbose,
    )
    return write_mrpt2_binary(filename, h1e, chunks, verbose=verbose, **kwargs)


# ---------------------------------------------------------------------------
# small demo
# ---------------------------------------------------------------------------


def _demo():
    import pyscf

    from pyscf_util.Integrals import integral_MRPT2_binary as this

    mol = pyscf.gto.Mole()
    mol.atom = "O 0 0 0; H 0 0.757 0.586; H 0 -0.757 0.586"
    mol.basis = "cc-pvdz"
    mol.symmetry = "C2v"
    mol.verbose = 0
    mol.build()
    mf = pyscf.scf.RHF(mol).run()
    mo = mf.mo_coeff
    nfzc, nact = 1, 6
    nvir = mol.nao - nfzc - nact
    tmpdir = tempfile.mkdtemp(prefix="mrpt2_bin_demo_")
    binfile = os.path.join(tmpdir, "dump.bin")
    stats = this.dump_mrpt2_binary(
        mol, mf, mo, nfzc, nact, nvir, binfile, 1e-10
    )
    data = this.load_mrpt2_binary(binfile)
    print(
        "nint=%d size=%.2fMB nint(check)=%d nmo=%d"
        % (
            stats["nint"],
            stats["bytes"] / 2.0**20,
            data["nint"],
            data["nmo"],
        )
    )
    return stats


if __name__ == "__main__":
    _demo()
