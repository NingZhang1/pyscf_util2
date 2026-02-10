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


def fcidump_size_consistency(
    # Mol1 Info #
    h1e_1,
    h2e_1,
    e_core_1,  #
    ncore_1,
    nact_1,
    nvirt_1,  #
    n_electron_1,  #
    # Mol2 Info #
    h1e_2,
    h2e_2,
    e_core_2,  #
    ncore_2,
    nact_2,
    nvirt_2,  #
    n_electron_2,
    _filename=None,
):
    # check #

    norb_1 = ncore_1 + nact_1 + nvirt_1
    norb_2 = ncore_2 + nact_2 + nvirt_2

    if len(h2e_1.shape) != 4:
        h2e_1 = pyscf.ao2mo.restore(1, h2e_1.copy(), norb_1)

    if len(h2e_2.shape) != 4:
        h2e_2 = pyscf.ao2mo.restore(1, h2e_2.copy(), norb_2)

    # build the full info #

    norb = norb_1 + norb_2
    ncore = ncore_1 + ncore_2
    nact = nact_1 + nact_2
    nvirt = nvirt_1 + nvirt_2
    n_electron = n_electron_1 + n_electron_2
    e_core = e_core_1 + e_core_2

    # build map #

    map_1 = []

    for i in range(ncore_1):
        map_1.append(i)

    for i in range(nact_1):
        map_1.append(ncore + i)

    for i in range(nvirt_1):
        map_1.append(ncore + nact + i)

    map_2 = []

    for i in range(ncore_2):
        map_2.append(ncore_1 + i)

    for i in range(nact_2):
        map_2.append(ncore + nact_1 + i)

    for i in range(nvirt_2):
        map_2.append(ncore + nact + nvirt_1 + i)

    # print(map_1)
    # print(map_2)

    # build h1e #

    h1e = np.zeros((norb, norb))

    for p in range(norb_1):
        for q in range(norb_1):
            h1e[map_1[p], map_1[q]] = h1e_1[p, q]

    for p in range(norb_2):
        for q in range(norb_2):
            h1e[map_2[p], map_2[q]] = h1e_2[p, q]

    # build h2e #

    h2e = np.zeros((norb, norb, norb, norb))

    for p in range(norb_1):
        for q in range(norb_1):
            for r in range(norb_1):
                for s in range(norb_1):
                    h2e[map_1[p], map_1[q], map_1[r], map_1[s]] = h2e_1[p, q, r, s]

    for p in range(norb_2):
        for q in range(norb_2):
            for r in range(norb_2):
                for s in range(norb_2):
                    h2e[map_2[p], map_2[q], map_2[r], map_2[s]] = h2e_2[p, q, r, s]

    orb_sym = [0 for _ in range(norb)]

    if _filename == None:
        return (
            h1e,
            h2e,
            e_core,
            norb,
            n_electron,
            orb_sym,
        )

    else:
        tools.fcidump.from_integrals(
            filename=_filename,
            h1e=h1e,
            h2e=h2e,
            nuc=e_core,
            nmo=norb,
            nelec=n_electron,  # Useless
            tol=1e-12,
            orbsym=orb_sym,
        )


def gfock_size_consistency(
    # Mol1 Info #
    gfock_1,
    ncore_1,
    nact_1,
    nvirt_1,  #
    # Mol2 Info #
    gfock_2,
    ncore_2,
    nact_2,
    nvirt_2,  #
):

    # check #

    norb_1 = ncore_1 + nact_1 + nvirt_1
    norb_2 = ncore_2 + nact_2 + nvirt_2

    # build the full info #

    norb = norb_1 + norb_2
    ncore = ncore_1 + ncore_2
    nact = nact_1 + nact_2
    nvirt = nvirt_1 + nvirt_2
    gfock = np.zeros((norb, norb))

    # build map #

    map_1 = []

    for i in range(ncore_1):
        map_1.append(i)

    for i in range(nact_1):
        map_1.append(ncore + i)

    for i in range(nvirt_1):
        map_1.append(ncore + nact + i)

    map_2 = []

    for i in range(ncore_2):
        map_2.append(ncore_1 + i)

    for i in range(nact_2):
        map_2.append(ncore + nact_1 + i)

    for i in range(nvirt_2):
        map_2.append(ncore + nact + nvirt_1 + i)

    # build gfock #

    for p in range(norb_1):
        for q in range(norb_1):
            gfock[map_1[p], map_1[q]] = gfock_1[p, q]

    for p in range(norb_2):
        for q in range(norb_2):
            gfock[map_2[p], map_2[q]] = gfock_2[p, q]

    return gfock
