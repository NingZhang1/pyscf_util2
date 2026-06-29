from pyscf import gto, scf, lib
import numpy
from functools import reduce
from pyscf import lib
from pyscf import tools
import pyscf

### some data ###

PARITY = {
    "s": 0,
    "p": 1,
    "d": 0,
    "f": 1,
    "g": 0,
    "h": 1,
    "i": 0,
}

CHARACTER_TABLE = {
    "D2h": {
        # E,C2x,C2y,C2z,i, sx,sy,sz,  \bar{E}
        "Ag": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        "B1g": [1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1],
        "B2g": [1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1],
        "B3g": [1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1],
        "Au": [1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1],
        "B1u": [1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1],
        "B2u": [1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1],
        "B3u": [1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1],
    }
}

Alpha_Beta_RepMat = {
    "D2h": [
        numpy.array([[1, 0], [0, 1]], dtype=numpy.complex128),  # E
        numpy.array([[0, -1.0j], [-1.0j, 0]], dtype=numpy.complex128),  # C2x
        numpy.array([[0, -1], [1, 0]], dtype=numpy.complex128),  # C2y
        numpy.array([[-1.0j, 0], [0, 1.0j]], dtype=numpy.complex128),  # C2z
        numpy.array([[1, 0], [0, 1]], dtype=numpy.complex128),  # i
        numpy.array([[0, -1.0j], [-1.0j, 0]], dtype=numpy.complex128),  # sx
        numpy.array([[0, -1], [1, 0]], dtype=numpy.complex128),  # sy
        numpy.array([[-1.0j, 0], [0, 1.0j]], dtype=numpy.complex128),  # sz
        numpy.array([[-1, 0], [0, -1]], dtype=numpy.complex128),  # \bar{E}
        numpy.array([[0, 1.0j], [1.0j, 0]], dtype=numpy.complex128),  # C2xE
        numpy.array([[0, 1], [-1, 0]], dtype=numpy.complex128),  # C2yE
        numpy.array([[1.0j, 0], [0, -1.0j]], dtype=numpy.complex128),  # C2zE
        numpy.array([[-1, 0], [0, -1]], dtype=numpy.complex128),  # iE
        numpy.array([[0, 1.0j], [1.0j, 0]], dtype=numpy.complex128),  # sxE
        numpy.array([[0, 1], [-1, 0]], dtype=numpy.complex128),  # syE
        numpy.array([[1.0j, 0], [0, -1.0j]], dtype=numpy.complex128),  # szE
    ]
}

IRREP_ID = {
    "D2h": {
        0: "Ag",
        1: "B1g",
        2: "B2g",
        3: "B3g",
        4: "Au",
        5: "B1u",
        6: "B2u",
        7: "B3u",
    },
}

### some utils ###


def get_rep_mat(symm_orb_id):
    norb = len(symm_orb_id)

    rep_mat = numpy.zeros((16, norb, norb), dtype=numpy.complex128)

    for i in range(norb):
        ID = IRREP_ID["D2h"][symm_orb_id[i]]
        character = CHARACTER_TABLE["D2h"][ID]
        for j in range(16):
            rep_mat[j, i, i] = character[j]

    return rep_mat


def get_alpha_beta_rep_mat():
    return numpy.array(Alpha_Beta_RepMat["D2h"])


def kron_prod(A, B):
    assert A.shape[0] == B.shape[0]
    C = []
    for i in range(A.shape[0]):
        C.append(numpy.kron(A[i], B[i]))
    return numpy.array(C)


### time reversal adaptation ###


def _apply_time_reversal_op(mol, mo_coeff, debug=False):
    """Calculate the time reversal operator in the basis of the mo_coeff

    Args:
        mol: a molecule object
        mo_coeff: the molecular orbital coefficients

    Kwargs:

    Returns:
        tr_act_packed: a list of [index, coefficient] for the time reversal operator

    """

    trmaps = mol.time_reversal_map()
    idxA = numpy.where(trmaps > 0)[0]
    idxB = trmaps[idxA] - 1
    n = trmaps.size
    idx2 = numpy.hstack((idxA, idxA + n, idxB, idxB + n))

    if debug:
        print("trmaps = ", trmaps)
        print("idxA   = ", idxA)
        print("idxB   = ", idxB)

    time_reversal_m = numpy.zeros((2 * n, 2 * n), dtype=numpy.int64)

    for irow, data in enumerate(trmaps):
        icol = data
        elmt = 1
        if data < 0:
            icol = -data - 1
            elmt = -1
        else:
            icol = data - 1
            elmt = 1
        time_reversal_m[irow, icol] = elmt

    time_reversal_m[n:, n:] = time_reversal_m[:n, :n]

    if debug:
        print("time_reversal_m = ", time_reversal_m)

    ovlp_4C = pyscf.scf.dhf.get_ovlp(mol)

    tr_act = reduce(
        numpy.dot, (mo_coeff.T.conj(), ovlp_4C, time_reversal_m, mo_coeff.conj())
    )
    tr_act_packed = []
    for i in range(tr_act.shape[0]):
        for j in range(tr_act.shape[1]):
            if abs(tr_act[i, j]) > 1e-6:
                if debug:
                    print("tr_act = ", i, j, tr_act[i, j])
                tr_act_packed.append([j, tr_act[i, j]])

    return tr_act_packed


def _time_reversal_symmetry_adapted(mol, mo_coeff, debug=False):
    """Adapt the molecular orbital coefficients to the time reversal symmetry

    Args:
        mol: a molecule object
        mo_coeff: the molecular orbital coefficients

    Kwargs:
        debug: whether to print the details

    Returns:

    """

    trmaps = mol.time_reversal_map()
    idxA = numpy.where(trmaps > 0)[0]
    idxB = trmaps[idxA] - 1
    n = trmaps.size
    idx2 = numpy.hstack((idxA, idxA + n, idxB, idxB + n))

    if debug:
        print("trmaps = ", trmaps)
        print("idxA   = ", idxA)
        print("idxB   = ", idxB)

    time_reversal_m = numpy.zeros((2 * n, 2 * n), dtype=numpy.int64)

    for irow, data in enumerate(trmaps):
        icol = data
        elmt = 1
        if data < 0:
            icol = -data - 1
            elmt = -1
        else:
            icol = data - 1
            elmt = 1
        time_reversal_m[irow, icol] = elmt

    time_reversal_m[n:, n:] = time_reversal_m[:n, :n]

    if debug:
        print("time_reversal_m = ", time_reversal_m)

    ovlp_4C = pyscf.scf.dhf.get_ovlp(mol)

    ######### the first step is to rotate the orb so that the TR is really to be addapted! #########

    Res = mo_coeff.copy()

    idxA_all = numpy.hstack((idxA, idxA + n))
    ovlp_A = ovlp_4C[idxA_all, :][:, idxA_all]

    for i in range(0, 2 * n, 2):
        mo_coeff_A = Res[idxA_all, i]
        norm_A = reduce(numpy.dot, (mo_coeff_A.T.conj(), ovlp_A, mo_coeff_A))
        if norm_A < 0.5:
            if debug:
                # print("norm_A = ", norm_A)
                # swap i and i+1
                print("swap %d and %d" % (i, i + 1))
            # Res[:, i] = -Res[:, i]
            tmp = Res[:, i].copy()
            Res[:, i] = Res[:, i + 1]
            Res[:, i + 1] = tmp
            # if debug:
            #     print(Res[idx2, i:i+2])

        ### real orbital ###

        mo_coeff_A = Res[:, i]
        real_A = mo_coeff_A.real
        norm_real = numpy.linalg.norm(real_A)
        imag_A = mo_coeff_A.imag
        norm_imag = numpy.linalg.norm(imag_A)
        if norm_imag > norm_real:
            print("times i at ", i)
            print("norm_imag = ", norm_imag)
            print("norm_real = ", norm_real)
            Res[:, i] = -1.0j * Res[:, i]

        mo_coeff_B = Res[:, i + 1]
        real_B = mo_coeff_B.real
        norm_real = numpy.linalg.norm(real_B)
        imag_B = mo_coeff_B.imag
        norm_imag = numpy.linalg.norm(imag_B)
        if norm_imag > norm_real:
            print("times i at ", i + 1)
            print("norm_imag = ", norm_imag)
            print("norm_real = ", norm_real)
            Res[:, i + 1] = -1.0j * Res[:, i + 1]

    ######### the second step is to reorder the orb so that the TR seems to be addapted! #########

    tr_act = reduce(numpy.dot, (Res.T.conj(), ovlp_4C, time_reversal_m, Res.conj()))
    tr_act_packed = []
    for i in range(tr_act.shape[0]):
        for j in range(tr_act.shape[1]):
            if abs(tr_act[i, j]) > 1e-6:
                if debug:
                    print("tr_act = ", i, j, tr_act[i, j])
                tr_act_packed.append([j, tr_act[i, j]])

    for i in range(0, 2 * n, 2):
        if tr_act_packed[i][0] != i + 1:
            print("Error in time reversal symmetry")
            exit(1)
        if tr_act_packed[i][1] < 0.0:
            if debug:
                print("plus -1 between %d and %d" % (i, i + 1))
            Res[:, i + 1] *= -1.0

    return Res


### Jz and d2h adaptation ###


def _atom_spinor_spatial_parity(mol, mo_coeff, debug=False):
    """Calculate the parity of the molecular orbitals

    Args:
        mol: a molecule object
        mo_coeff: the molecular orbital coefficients

    Kwargs:
        debug: whether to print the details

    Returns:
        Res: the parity of the molecular orbitals, a list of 0 or 1, 0 for even, 1 for odd

    """

    labels = mol.spinor_labels()
    if debug:
        print(labels)

    n = mol.nao_2c()
    parity = []
    for i in range(n):
        for key in PARITY.keys():
            if key in labels[i]:
                parity.append(PARITY[key])
                break

    if debug:
        print("parity = ", parity)

    # add 4C's partity

    def _reverse_parity(parity):
        if parity == 0:
            return 1
        else:
            return 0

    # parity_small = [_reverse_parity(x) for x in parity]

    parity.extend(parity)

    # extract parity of the mo_coeff

    # the order of AO is first large component then small component

    ovlp_4C = pyscf.scf.dhf.get_ovlp(mol)

    pairty_even = [id for id, x in enumerate(parity) if x == 0]
    parity_odd = [id for id, x in enumerate(parity) if x == 1]

    ovlp_even = ovlp_4C[pairty_even, :][:, pairty_even]
    ovlp_odd = ovlp_4C[parity_odd, :][:, parity_odd]

    if debug:
        ovlp_cross = ovlp_4C[pairty_even, :][:, parity_odd]  # should all be zero
        print("ovlp_cross should be zero")
        print(ovlp_cross.shape)
        # print(ovlp_cross)
        print(numpy.allclose(ovlp_cross, numpy.zeros_like(ovlp_cross)))
        # print("pairty_even = ", pairty_even)

    # loop each orb and find their parity

    Res = []
    for i in range(mo_coeff.shape[1]):
        mo_coeff_A = mo_coeff[:, i].reshape(-1, 1)
        # print("mo_coeff_A.shape = ", mo_coeff_A.shape)
        # print("mo_coeff_A[pairty_even] = ", mo_coeff_A[pairty_even].shape)
        norm_even = reduce(
            numpy.dot,
            (mo_coeff_A[pairty_even].T.conj(), ovlp_even, mo_coeff_A[pairty_even]),
        )
        norm_odd = reduce(
            numpy.dot,
            (mo_coeff_A[parity_odd].T.conj(), ovlp_odd, mo_coeff_A[parity_odd]),
        )
        if norm_even < norm_odd:
            Res.append(1)
        else:
            Res.append(0)
        if debug:
            print("norm_even = ", norm_even, " norm_odd = ", norm_odd)
            print("Res = ", Res[-1])
            assert numpy.allclose(norm_even, 0.0) or numpy.allclose(norm_odd, 0.0)
    return Res


def _get_Jz_AO(mol, debug=False):

    labels = mol.spinor_labels()
    if debug:
        print(labels)
    n = mol.nao_2c()
    Jz = numpy.zeros((2 * n, 2 * n), dtype=numpy.complex128)
    # parity = []
    for i in range(n):
        sz_str = labels[i].split(",")[1]
        for j in range(len(sz_str)):
            if sz_str[j] == " ":
                sz_str = sz_str[:j]
                break
        a, b = sz_str.split("/")
        Jz[i, i] = float(int(a) / int(b))
        Jz[i + n, i + n] = float(int(a) / int(b))

    return Jz


def _atom_Jz_adapted(mol, mo_coeff, mo_energy, debug=False):
    """Adapt the molecular orbital coefficients to the Jz symmetry

    Args:
        mol: a molecule object, with only one atom
        mo_coeff: the molecular orbital coefficients
        mo_energy: the molecular orbital energies

    Kwargs:
        debug: whether to print the details

    Returns:

    """

    assert mol.natm == 1

    labels = mol.spinor_labels()
    if debug:
        print(labels)
    n = mol.nao_2c()
    Jz = numpy.zeros((2 * n, 2 * n), dtype=numpy.complex128)

    for i in range(n):
        sz_str = labels[i].split(",")[1]
        for j in range(len(sz_str)):
            if sz_str[j] == " ":
                sz_str = sz_str[:j]
                break
        a, b = sz_str.split("/")
        Jz[i, i] = float(int(a) / int(b))
        Jz[i + n, i + n] = float(int(a) / int(b))

    if debug:
        print("Jz     = ", Jz)

    # the order of AO is first large component then small component

    ovlp_4C = pyscf.scf.dhf.get_ovlp(mol)

    loc = 0

    permute = []

    Res = mo_coeff.copy()

    while True:

        ene_now = mo_energy[loc]
        loc_end = None
        for i in range(loc, mo_coeff.shape[1]):
            # if abs(mo_energy[i] - ene_now) > 1e-5:
            # if abs(mo_energy[i] - ene_now) > 1e-7:
            if abs(mo_energy[i] - ene_now) > 1e-8:  # 判断哪些轨道属于同一个能级
                loc_end = i
                break

        if debug:
            print("loc = ", loc, " loc_end = ", loc_end, " ene_now = ", ene_now)
            print("mo_energy[loc:loc_end] = ", mo_energy[loc:loc_end])

        if loc_end is None:
            loc_end = mo_coeff.shape[1]

        norb = loc_end - loc

        mo_coeff_tmp = mo_coeff[:, loc:loc_end]
        Jz_Tmp = reduce(numpy.dot, (mo_coeff_tmp.T.conj(), ovlp_4C, Jz, mo_coeff_tmp))
        e, h = numpy.linalg.eigh(Jz_Tmp)

        if debug:
            print("Jz_Tmp = ", Jz_Tmp)
            print("e = ", e)

        # rotate the orbitals

        Res[:, loc:loc_end] = numpy.dot(mo_coeff_tmp, h)

        for i in range(norb // 2):
            permute.append(loc + norb - 1 - i)
            permute.append(loc + i)

        if i % 2 == 1:
            if debug:
                print("--------------------")

        loc = loc_end

        if debug:
            print("loc = ", loc)
            print("mf.mo_coeff.shape[1] = ", Res.shape[1])

        if loc == mo_coeff.shape[1]:
            break

    if debug:
        print("permute = ", permute)

    Res = Res[:, permute]

    if debug:
        ovlp_mo = reduce(numpy.dot, (Res.T.conj(), ovlp_4C, Res))
        print("ovlp_mo = ", numpy.diag(ovlp_mo))

        Jz_mo = reduce(numpy.dot, (Res.T.conj(), ovlp_4C, Jz, Res))
        print("Jz_mo = ", numpy.diag(Jz_mo)[n:])

    return Res


def _atm_spinor_2_d2h_adapted_spinor(mol):

    TRANS_MAP = {
        # J = p + 1/2, p odd, Jz = q + 0.5 q odd, in pyscf the order is |J-Jz>, |J Jz>
        0: numpy.array([[-1, 0], [0, 1]]),
        # J = p + 1/2, p odd, Jz = q + 0.5 q even, in pyscf the order is |J-Jz>, |J Jz>
        1: numpy.array([[0, -1], [1, 0]]),
        # J = p + 1/2, p even, Jz = q + 0.5 q odd, in pyscf the order is |J-Jz>, |J Jz>
        2: numpy.array([[1, 0], [0, 1]]),
        # J = p + 1/2, p even, Jz = q + 0.5 q even, in pyscf the order is |J-Jz>, |J Jz>
        3: numpy.array([[0, 1], [1, 0]]),
    }

    LOC_MAP = {
        1: [3],
        3: [0, 1],
        5: [3, 2, 3],
        7: [0, 1, 0, 1],
        9: [3, 2, 3, 2, 3],
        11: [0, 1, 0, 1, 0, 1],
        13: [3, 2, 3, 2, 3, 2, 3],
        15: [0, 1, 0, 1, 0, 1, 0, 1],
    }

    n2c = mol.nao_2c()

    Res = numpy.zeros((n2c, n2c), dtype=numpy.complex128)

    indxA = []

    loc = 0
    loc_col = 0
    for ib in range(mol.nbas):
        l = mol.bas_angular(ib)
        kappa = mol.bas_kappa(ib)
        print("basis %3d l = %3d kappa = %3d" % (ib, l, kappa))
        assert kappa == 0
        nctr = mol.bas_nctr(ib)
        print("nctr = ", nctr)
        if l == 0:
            for _ in range(nctr):
                # p = 0 and q = 0
                trans_map = TRANS_MAP[3]
                Res[loc, loc_col] = trans_map[0, 0]
                Res[loc, loc_col + 1] = trans_map[0, 1]
                Res[loc + 1, loc_col] = trans_map[1, 0]
                Res[loc + 1, loc_col + 1] = trans_map[1, 1]
                indxA.append(loc_col)
                loc += 2
                loc_col += 2
        elif l == 1:
            for _ in range(nctr):
                # J 1/2 Jz 1/2
                # p = 0 and q = 0
                trans_map = TRANS_MAP[3]
                Res[loc, loc_col] = trans_map[0, 0]
                Res[loc, loc_col + 1] = trans_map[0, 1]
                Res[loc + 1, loc_col] = trans_map[1, 0]
                Res[loc + 1, loc_col + 1] = trans_map[1, 1]
                indxA.append(loc_col)
                loc_col += 2
                loc += 2
                # # J 3/2 Jz 3/2
                # # p = 1 and q = 1
                # trans_map = TRANS_MAP[0]
                # Res[loc, loc] = trans_map[0, 0]
                # Res[loc, loc+4-1] = trans_map[0, 1]
                # Res[loc+4-1, loc] = trans_map[1, 0]
                # Res[loc+4-1, loc+4-1] = trans_map[1, 1]
                # # J 3/2 Jz 1/2
                # # p = 1 and q = 0
                # trans_map = TRANS_MAP[1]
                # Res[loc+1, loc+1] = trans_map[0, 0]
                # Res[loc+1, loc+2] = trans_map[0, 1]
                # Res[loc+2, loc+1] = trans_map[1, 0]
                # Res[loc+2, loc+2] = trans_map[1, 1]

                for loc_tmp, info in enumerate(LOC_MAP[3]):
                    trans_map = TRANS_MAP[info]
                    Res[loc + loc_tmp, loc_col] = trans_map[0, 0]
                    Res[loc + loc_tmp, loc_col + 1] = trans_map[0, 1]
                    Res[loc + 3 - loc_tmp, loc_col] = trans_map[1, 0]
                    Res[loc + 3 - loc_tmp, loc_col + 1] = trans_map[1, 1]
                    indxA.append(loc_col)
                    loc_col += 2
                loc += 4
        elif l == 2:
            for _ in range(nctr):
                # J 3/2
                for loc_tmp, info in enumerate(LOC_MAP[3]):
                    trans_map = TRANS_MAP[info]
                    Res[loc + loc_tmp, loc_col] = trans_map[0, 0]
                    Res[loc + loc_tmp, loc_col + 1] = trans_map[0, 1]
                    Res[loc + 3 - loc_tmp, loc_col] = trans_map[1, 0]
                    Res[loc + 3 - loc_tmp, loc_col + 1] = trans_map[1, 1]
                    indxA.append(loc_col)
                    loc_col += 2
                loc += 4
                # J 5/2
                for loc_tmp, info in enumerate(LOC_MAP[5]):
                    trans_map = TRANS_MAP[info]
                    Res[loc + loc_tmp, loc_col] = trans_map[0, 0]
                    Res[loc + loc_tmp, loc_col + 1] = trans_map[0, 1]
                    Res[loc + 5 - loc_tmp, loc_col] = trans_map[1, 0]
                    Res[loc + 5 - loc_tmp, loc_col + 1] = trans_map[1, 1]
                    indxA.append(loc_col)
                    loc_col += 2
                loc += 6
        elif l == 3:
            for _ in range(nctr):
                # J 5/2
                for loc_tmp, info in enumerate(LOC_MAP[5]):
                    trans_map = TRANS_MAP[info]
                    Res[loc + loc_tmp, loc_col] = trans_map[0, 0]
                    Res[loc + loc_tmp, loc_col + 1] = trans_map[0, 1]
                    Res[loc + 5 - loc_tmp, loc_col] = trans_map[1, 0]
                    Res[loc + 5 - loc_tmp, loc_col + 1] = trans_map[1, 1]
                    indxA.append(loc_col)
                    loc_col += 2
                loc += 6
                # J 7/2
                for loc_tmp, info in enumerate(LOC_MAP[7]):
                    trans_map = TRANS_MAP[info]
                    Res[loc + loc_tmp, loc_col] = trans_map[0, 0]
                    Res[loc + loc_tmp, loc_col + 1] = trans_map[0, 1]
                    Res[loc + 7 - loc_tmp, loc_col] = trans_map[1, 0]
                    Res[loc + 7 - loc_tmp, loc_col + 1] = trans_map[1, 1]
                    indxA.append(loc_col)
                    loc_col += 2
                loc += 8
        elif l == 4:
            for _ in range(nctr):
                # J 7/2
                for loc_tmp, info in enumerate(LOC_MAP[7]):
                    trans_map = TRANS_MAP[info]
                    Res[loc + loc_tmp, loc_col] = trans_map[0, 0]
                    Res[loc + loc_tmp, loc_col + 1] = trans_map[0, 1]
                    Res[loc + 7 - loc_tmp, loc_col] = trans_map[1, 0]
                    Res[loc + 7 - loc_tmp, loc_col + 1] = trans_map[1, 1]
                    indxA.append(loc_col)
                    loc_col += 2
                loc += 8
                # J 9/2
                for loc_tmp, info in enumerate(LOC_MAP[9]):
                    trans_map = TRANS_MAP[info]
                    Res[loc + loc_tmp, loc_col] = trans_map[0, 0]
                    Res[loc + loc_tmp, loc_col + 1] = trans_map[0, 1]
                    Res[loc + 9 - loc_tmp, loc_col] = trans_map[1, 0]
                    Res[loc + 9 - loc_tmp, loc_col + 1] = trans_map[1, 1]
                    indxA.append(loc_col)
                    loc_col += 2
                loc += 10
        elif l == 5:
            for _ in range(nctr):
                # J 9/2
                for loc_tmp, info in enumerate(LOC_MAP[9]):
                    trans_map = TRANS_MAP[info]
                    Res[loc + loc_tmp, loc_col] = trans_map[0, 0]
                    Res[loc + loc_tmp, loc_col + 1] = trans_map[0, 1]
                    Res[loc + 9 - loc_tmp, loc_col] = trans_map[1, 0]
                    Res[loc + 9 - loc_tmp, loc_col + 1] = trans_map[1, 1]
                    indxA.append(loc_col)
                    loc_col += 2
                loc += 10
                # J 11/2
                for loc_tmp, info in enumerate(LOC_MAP[11]):
                    trans_map = TRANS_MAP[info]
                    Res[loc + loc_tmp, loc_col] = trans_map[0, 0]
                    Res[loc + loc_tmp, loc_col + 1] = trans_map[0, 1]
                    Res[loc + 11 - loc_tmp, loc_col] = trans_map[1, 0]
                    Res[loc + 11 - loc_tmp, loc_col + 1] = trans_map[1, 1]
                    indxA.append(loc_col)
                    loc_col += 2
                loc += 12
        elif l == 6:
            for _ in range(nctr):
                # J 11/2
                for loc_tmp, info in enumerate(LOC_MAP[11]):
                    trans_map = TRANS_MAP[info]
                    Res[loc + loc_tmp, loc_col] = trans_map[0, 0]
                    Res[loc + loc_tmp, loc_col + 1] = trans_map[0, 1]
                    Res[loc + 11 - loc_tmp, loc_col] = trans_map[1, 0]
                    Res[loc + 11 - loc_tmp, loc_col + 1] = trans_map[1, 1]
                    indxA.append(loc_col)
                    loc_col += 2
                loc += 12
                # J 13/2
                for loc_tmp, info in enumerate(LOC_MAP[13]):
                    trans_map = TRANS_MAP[info]
                    Res[loc + loc_tmp, loc_col] = trans_map[0, 0]
                    Res[loc + loc_tmp, loc_col + 1] = trans_map[0, 1]
                    Res[loc + 13 - loc_tmp, loc_col] = trans_map[1, 0]
                    Res[loc + 13 - loc_tmp, loc_col + 1] = trans_map[1, 1]
                    indxA.append(loc_col)
                    loc_col += 2
                loc += 14
        else:
            raise ValueError("l = %3d is not supported" % l)

    return Res, indxA


### time rversal symm adapted ###


def time_reversal_symm_adapted(mol, mo_coeff):

    n2c = mol.nao_2c()
    if mo_coeff.shape[1] == n2c:
        mo_pes = mo_coeff.copy()
    else:
        mo_pes = mo_coeff[:, n2c:]

    tr_mat = _apply_time_reversal_op(mol, mo_pes, debug=True)

    for i in range(0, mo_pes.shape[1], 2):
        assert tr_mat[i][0] == i + 1
        if abs(tr_mat[i][1].real + 1) < 1e-4:
            if True:
                print("swap the sign of the orbitals ", i + 1)
            mo_pes[:, i + 1] *= -1.0

    return mo_pes


### the only working function ###


def atm_d2h_symmetry_adapt_mo_coeff(mol, mo_coeff, debug=False):

    spinor_2_adapted, indA = _atm_spinor_2_d2h_adapted_spinor(mol)
    n2c = mol.nao_2c()
    indB = [i for i in range(n2c) if i not in indA]
    if debug:
        print("indA = ", indA)
        print("indB = ", indB)

    if mo_coeff.shape[1] == n2c:
        mo_pes = mo_coeff.copy()
    else:
        mo_pes = mo_coeff[:, n2c:]
    mo_pes_large_comp = mo_pes[:n2c, :]

    mo_coeff_over_adapted_spinor = numpy.dot(
        spinor_2_adapted.conj().T, mo_pes_large_comp
    )

    if debug:
        for j in range(mo_coeff_over_adapted_spinor.shape[1]):
            for i in range(mo_coeff_over_adapted_spinor.shape[0]):
                if abs(mo_coeff_over_adapted_spinor[i, j]) > 1e-8:
                    print(
                        "%3d %3d %20.10f %20.10f"
                        % (
                            i,
                            j,
                            mo_coeff_over_adapted_spinor[i, j].real,
                            mo_coeff_over_adapted_spinor[i, j].imag,
                        )
                    )

    ##################### first step make sure that all the orbitals are real #####################

    for i in range(mo_coeff_over_adapted_spinor.shape[1]):
        coeff_tmp = mo_coeff_over_adapted_spinor[:, i]
        norm_real = numpy.linalg.norm(coeff_tmp.real)
        norm_imag = numpy.linalg.norm(coeff_tmp.imag)
        if debug:
            print("norm_real = ", norm_real, "norm_imag = ", norm_imag)
        if norm_real < norm_imag:
            mo_pes[:, i] *= 1.0j

    mo_pes_large_comp = mo_pes[:n2c, :]
    mo_coeff_over_adapted_spinor = numpy.dot(
        spinor_2_adapted.conj().T, mo_pes_large_comp
    ).real

    ##################### second step make sure that all the orbitals are real #####################

    # time reversal adapted

    # first switch the comp

    for i in range(0, mo_coeff_over_adapted_spinor.shape[1], 2):
        coeff_tmp1 = mo_coeff_over_adapted_spinor[:, i]

        norm_A = numpy.linalg.norm(coeff_tmp1[indA])
        norm_B = numpy.linalg.norm(coeff_tmp1[indB])
        if debug:
            print("norm_A = ", norm_A, "norm_B = ", norm_B)
        if norm_A < norm_B:
            if debug:
                print("swap the order of the orbitals ", i, i + 1)
            # swap the order of the orbitals
            import copy

            coeff_tmp1 = copy.deepcopy(mo_coeff_over_adapted_spinor[:, i])
            mo_coeff_over_adapted_spinor[:, i] = mo_coeff_over_adapted_spinor[:, i + 1]
            mo_coeff_over_adapted_spinor[:, i + 1] = coeff_tmp1
            # swap mo_pes
            tmp = copy.deepcopy(mo_pes[:, i])
            mo_pes[:, i] = mo_pes[:, i + 1]
            mo_pes[:, i + 1] = tmp

    # second the sign problem

    # from Util_Rela4C import _apply_time_reversal_op

    tr_mat = _apply_time_reversal_op(mol, mo_pes, debug=debug)

    # print("mo_pes.shape = ", mo_pes.shape)

    for i in range(0, mo_pes.shape[1], 2):
        assert tr_mat[i][0] == i + 1
        if abs(tr_mat[i][1].real + 1) < 1e-4:
            if debug:
                print("swap the sign of the orbitals ", i + 1)
            mo_pes[:, i + 1] *= -1.0

    # if debug:
    #     _apply_time_reversal_op(mol, mo_pes, debug=debug)

    # generate new mo_coeff

    if debug:
        # check
        mo_pes_large_comp = mo_pes[:n2c, :]
        mo_coeff_over_adapted_spinor = numpy.dot(
            spinor_2_adapted.conj().T, mo_pes_large_comp
        )

        for j in range(mo_coeff_over_adapted_spinor.shape[1]):
            for i in range(mo_coeff_over_adapted_spinor.shape[0]):
                if abs(mo_coeff_over_adapted_spinor[i, j]) > 1e-8:
                    print(
                        "%3d %3d %20.10f %20.10f"
                        % (
                            i,
                            j,
                            mo_coeff_over_adapted_spinor[i, j].real,
                            mo_coeff_over_adapted_spinor[i, j].imag,
                        )
                    )

    return mo_pes


if __name__ == "__main__":
    # mol = gto.M(atom='H 0 0 0; H 0 0 1; O 0 1 0', basis='sto-3g', verbose=5)
    mol = gto.M(
        atom="F 0 0 0", basis="cc-pvdz", verbose=5, charge=-1, spin=0, symmetry="d2h"
    )
    mol.build()
    mf = scf.dhf.RDHF(mol)
    mf.conv_tol = 1e-12
    # mf.with_breit = True
    mf.with_gaunt = True
    mf.kernel()

    mo_coeff_backup = mf.mo_coeff.copy()

    # the benchmark

    # _, _ = FCIDUMP_Rela4C(
    #     mol, mf, with_breit=True, filename="FCIDUMP_4C_Breit", mode="original", debug=False)
    # _, _ = FCIDUMP_Rela4C(
    #     mol, mf, with_breit=True, filename="FCIDUMP_4C_Breit2", mode="incore", debug=False)

    fock = mf.get_fock()

    print("mo_ene = ", mf.mo_energy)

    ovlp_4C = pyscf.scf.dhf.get_ovlp(mol)

    e, mo_coeff = mf._eigh(fock, ovlp_4C)

    print("mo_ene = ", e)

    print(numpy.allclose(mf.mo_energy, e))
    print(numpy.allclose(mf.mo_coeff, mo_coeff))

    mf.mo_coeff = mo_coeff

    # test the symmetry adapted basis

    spinor_2_adapted, indA = _atm_spinor_2_d2h_adapted_spinor(mol)

    print("indA = ", indA)

    print(mol.symm_orb)
    print(mol.irrep_id)

    sym_orb = numpy.concatenate([mol.symm_orb[0], mol.symm_orb[1]], axis=1)
    for i in range(2, len(mol.symm_orb)):
        sym_orb = numpy.concatenate([sym_orb, mol.symm_orb[i]], axis=1)

    from pyscf import symm

    orbsym_ID = symm.label_orb_symm(mol, mol.irrep_id, mol.symm_orb, numpy.eye(mol.nao))

    # print(get_rep_mat(orbsym_ID))
    rep_mat = get_rep_mat(orbsym_ID)
    for i in range(16):
        print(numpy.diag(rep_mat[i]))

    ca, cb = pyscf.symm.sph.sph2spinor_coeff(mol)
    A = numpy.concatenate([ca, cb], axis=0)
    print(A.shape)
    print(A[:, -1])

    # print non-zero elements

    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            if abs(A[i, j]) > 1e-8:
                print(i, j, A[i, j])

    A_adapted = numpy.dot(A, spinor_2_adapted)

    alpha_beta_repmat = get_alpha_beta_rep_mat()
    print(rep_mat.shape)
    print(alpha_beta_repmat.shape)

    kron_prod_mat = kron_prod(alpha_beta_repmat, rep_mat)

    # calculate the rep matrix over symmetrized spinor basis

    ovlp_2C = numpy.zeros((mol.nao * 2, mol.nao * 2), dtype=numpy.complex128)
    ovlp_2C[: mol.nao, : mol.nao] = mol.intor("int1e_ovlp_sph")
    ovlp_2C[mol.nao :, mol.nao :] = mol.intor("int1e_ovlp_sph")

    # check spinor ovlp

    ovlp_spinor = reduce(numpy.dot, (A.T.conj(), ovlp_2C, A))
    ovlp_spinor_bench = mol.intor_symmetric("int1e_ovlp_spinor")
    print("ovlp all the same ? ", numpy.allclose(ovlp_spinor, ovlp_spinor_bench))

    for i in range(16):

        continue

        rep_now = reduce(numpy.dot, (A_adapted.T.conj(), kron_prod_mat[i], A_adapted))

        print("Rep Mat for op ", i)

        # the correct way to print the rep mat

        loc = 0

        for ib in range(mol.nbas):
            l = mol.bas_angular(ib)
            kappa = mol.bas_kappa(ib)
            print("basis %3d l = %3d kappa = %3d" % (ib, l, kappa))
            assert kappa == 0
            nctr = mol.bas_nctr(ib)
            print("nctr = ", nctr)
            if l == 0:
                for _ in range(nctr):
                    print(rep_now[loc : loc + 2, :][:, loc : loc + 2])
                    loc += 2
            elif l == 1:
                for _ in range(nctr):
                    print(rep_now[loc : loc + 2, :][:, loc : loc + 2])
                    loc += 2
                    print(rep_now[loc : loc + 2, :][:, loc : loc + 2])
                    print(rep_now[loc + 2 : loc + 4, :][:, loc + 2 : loc + 4])
                    loc += 4
            elif l == 2:
                for _ in range(nctr):
                    print(rep_now[loc : loc + 2, :][:, loc : loc + 2])
                    print(rep_now[loc + 2 : loc + 4, :][:, loc + 2 : loc + 4])
                    loc += 4
                    print(rep_now[loc : loc + 2, :][:, loc : loc + 2])
                    print(rep_now[loc + 2 : loc + 4, :][:, loc + 2 : loc + 4])
                    print(rep_now[loc + 4 : loc + 6, :][:, loc + 4 : loc + 6])
                    loc += 6

    mf.mo_coeff = _atom_Jz_adapted(mol, mf.mo_coeff, mf.mo_energy, True)

    mo_pes = atm_d2h_symmetry_adapt_mo_coeff(mol, mf.mo_coeff, True)

    for j in range(mo_pes.shape[1]):
        for i in range(mo_pes.shape[0]):
            if abs(mo_pes[i, j]) > 1e-8:
                print(
                    "%3d %3d %20.10f %20.10f"
                    % (i, j, mo_pes[i, j].real, mo_pes[i, j].imag)
                )
                assert abs(mo_pes[i, j].imag) < 1e-8
