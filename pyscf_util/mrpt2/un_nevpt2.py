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
    Mol, my_mcscf, mo_coeff, orb_ene, ncore, nact, nvirt, _filename="FCIDUMP_Dyall"
):

    norb = ncore + nact + nvirt

    int2e_full = pyscf.ao2mo.full(
        eri_or_mol=Mol, mo_coeff=mo_coeff[:, ncore : ncore + nact], compact=True
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
        my_mcscf, mo_coeff=mo_coeff[:, : ncore + nact], ncas=nact, ncore=ncore
    )

    h1e = np.zeros((norb, norb))
    h1e[ncore : ncore + nact, ncore : ncore + nact] = int1e_res

    # set orb_ene #

    int1e_res = h1e
    for i in range(ncore):
        int1e_res[i, i] = orb_ene[i]
    for i in range(ncore + nact, norb):
        int1e_res[i, i] = orb_ene[i]

    # shift energy_core #

    for i in range(ncore):
        energy_core -= 2 * orb_ene[i]

    # get core #

    if _filename == None:
        return (
            int1e_res,
            int2e_res,
            energy_core,
            norb,
            Mol.nelectron,
            OrbSymID,
        )
    else:
        tools.fcidump.from_integrals(
            filename=_filename,
            h1e=int1e_res,
            h2e=int2e_res,
            nuc=energy_core,
            nmo=norb,
            nelec=Mol.nelectron,  # Useless
            tol=1e-10,
            orbsym=OrbSymID,
        )
