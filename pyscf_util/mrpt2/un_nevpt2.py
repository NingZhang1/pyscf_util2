# coding=UTF-8

import pyscf
from pyscf import tools
from pyscf import symm
from pyscf.tools import fcidump
import pyscf.mcscf
from pyscf_util.File import file_rdm
import numpy as np
from pyscf import tools
from pyscf.ao2mo import incore
from functools import reduce


def fcidump_Dyall(
    Mol,
    my_mcscf,
    mo_coeff,
    orb_ene,
    ncore,
    nact,
    nvirt,
    nfzc=0,
    _filename="FCIDUMP_Dyall",
):

    norb = ncore + nact + nvirt

    int2e_full = pyscf.ao2mo.full(
        eri_or_mol=Mol,
        mo_coeff=mo_coeff[:, nfzc + ncore : nfzc + ncore + nact],
        compact=True,
    )  # incore anyway since the size of active space cannot be too large!
    int2e_full = pyscf.ao2mo.restore(1, int2e_full.copy(), nact)

    int2e_res = np.zeros((norb, norb, norb, norb))
    int2e_res[
        ncore : ncore + nact,
        ncore : ncore + nact,
        ncore : ncore + nact,
        ncore : ncore + nact,
    ] = int2e_full

    # scf = my_mcscf._scf
    # h1e = reduce(np.dot, (mo_coeff.T, scf.get_hcore(), mo_coeff))

    # get orbsym

    OrbSym = pyscf.symm.label_orb_symm(Mol, Mol.irrep_name, Mol.symm_orb, mo_coeff)
    OrbSymID = [pyscf.symm.irrep_name2id(Mol.groupname, x) for x in OrbSym]

    int1e_res, energy_core = pyscf.mcscf.casci.h1e_for_cas(
        my_mcscf,
        mo_coeff=mo_coeff[:, : nfzc + ncore + nact],
        ncas=nact,
        ncore=nfzc + ncore,
    )

    h1e = np.zeros((norb, norb))

    # print(ncore, ncore + nact, nfzc + ncore, nfzc + ncore + nact)

    h1e[ncore : ncore + nact, ncore : ncore + nact] = int1e_res

    # set orb_ene #

    int1e_res = h1e
    for i in range(nfzc, nfzc + ncore):
        int1e_res[i - nfzc, i - nfzc] = orb_ene[i]
    for i in range(nfzc + ncore + nact, nfzc + norb):
        int1e_res[i - nfzc, i - nfzc] = orb_ene[i]

    # shift energy_core #

    for i in range(nfzc, nfzc + ncore):
        energy_core -= 2 * orb_ene[i]

    # for i in range(nfzc):
    #     energy_core += 2 * orb_ene[i]

    # get core #

    if _filename == None:
        return (
            int1e_res,
            int2e_res,
            energy_core,
            norb,
            Mol.nelectron,
            OrbSymID[nfzc : nfzc + norb],
        )
    else:
        tools.fcidump.from_integrals(
            filename=_filename,
            h1e=int1e_res,
            h2e=int2e_res,
            nuc=energy_core,
            nmo=norb,
            nelec=Mol.nelectron - 2 * nfzc,  # Useless
            tol=1e-10,
            orbsym=OrbSymID[nfzc : nfzc + norb],
        )


def fcidump_gFock(
    Mol,
    # my_mcscf,
    mo_coeff,
    gfock,
    # orb_ene,
    # ncore,
    # nact,
    # nvirt,
    nfzc=0,
    _filename="FCIDUMP_gFock",
):

    # norb = ncore + nact + nvirt

    norb = mo_coeff.shape[1] - nfzc

    # int2e_full = pyscf.ao2mo.full(
    #     eri_or_mol=Mol,
    #     mo_coeff=mo_coeff[:, nfzc + ncore : nfzc + ncore + nact],
    #     compact=True,
    # )  # incore anyway since the size of active space cannot be too large!
    # int2e_full = pyscf.ao2mo.restore(1, int2e_full.copy(), nact)

    # int2e_res = np.zeros((norb, norb, norb, norb))
    # int2e_res[
    #     ncore : ncore + nact,
    #     ncore : ncore + nact,
    #     ncore : ncore + nact,
    #     ncore : ncore + nact,
    # ] = int2e_full

    int2e_res = np.zeros((norb, norb, norb, norb))

    # scf = my_mcscf._scf
    # h1e = reduce(np.dot, (mo_coeff.T, scf.get_hcore(), mo_coeff))

    # get orbsym

    OrbSym = pyscf.symm.label_orb_symm(Mol, Mol.irrep_name, Mol.symm_orb, mo_coeff)
    OrbSymID = [pyscf.symm.irrep_name2id(Mol.groupname, x) for x in OrbSym]

    # int1e_res, energy_core = pyscf.mcscf.casci.h1e_for_cas(
    #     my_mcscf,
    #     mo_coeff=mo_coeff[:, : nfzc + ncore + nact],
    #     ncas=nact,
    #     ncore=nfzc + ncore,
    # )

    # h1e = np.zeros((norb, norb))

    # print(ncore, ncore + nact, nfzc + ncore, nfzc + ncore + nact)

    # h1e[ncore : ncore + nact, ncore : ncore + nact] = int1e_res

    h1e = gfock[nfzc:, nfzc:]

    # set orb_ene #

    int1e_res = h1e
    # for i in range(nfzc, nfzc + ncore):
    #     int1e_res[i - nfzc, i - nfzc] = orb_ene[i]
    # for i in range(nfzc + ncore + nact, nfzc + norb):
    #     int1e_res[i - nfzc, i - nfzc] = orb_ene[i]

    # shift energy_core #

    energy_core = 0.0
    for i in range(0, nfzc):
        energy_core += 2 * gfock[i, i]

    # for i in range(nfzc):
    #     energy_core += 2 * orb_ene[i]

    # get core #

    if _filename == None:
        return (
            int1e_res,
            int2e_res,
            energy_core,
            norb,
            Mol.nelectron,
            OrbSymID[nfzc:],
        )
    else:
        tools.fcidump.from_integrals(
            filename=_filename,
            h1e=int1e_res,
            h2e=int2e_res,
            nuc=energy_core,
            nmo=norb,
            nelec=Mol.nelectron - 2 * nfzc,  # Useless
            tol=1e-10,
            orbsym=OrbSymID[nfzc:],
        )
