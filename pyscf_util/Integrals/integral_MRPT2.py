import pyscf
import numpy
from functools import reduce
from pyscf import tools
import tempfile
import h5py
from pyscf.ao2mo import outcore
from pyscf.mcscf.casci import get_fock
from pyscf_util.misc.misc import _combine4, _combine2
from pyscf_util.File import file_rdm

############## build generalized Fock operator ##############

from functools import reduce
import numpy


def get_generalized_fock(
    mc, mo_coeff, rdm1, with_g3=False, only_g3=False, _opt_impl=True
):
    """
    casdm1 (ndarray): 1-particle density matrix in active space. Without
            input casdm1, the density matrix is computed with the input ci
            coefficients/object. If neither ci nor casdm1 were given, density
            matrix is computed by :func:`mc.fcisolver.make_rdm1` method. For
            state-average CASCI/CASCF calculation, this results in the
            effective Fock matrix based on the state-average density matrix.
            To obtain the effective Fock matrix for one particular state, you
            can assign the density matrix of that state to the kwarg casdm1.
    """

    assert rdm1 is not None
    assert mo_coeff is not None

    if only_g3 and not with_g3:
        nmo = mo_coeff.shape[1]
        return numpy.zeros((nmo, nmo))

    if not only_g3:
        fock_ao = get_fock(mc, mo_coeff=mo_coeff, casdm1=rdm1)

    if not with_g3:
        return reduce(numpy.dot, (mo_coeff.conj().T, fock_ao, mo_coeff))
    else:
        d = 2.0 * numpy.identity(rdm1.shape[0]) - rdm1
        # print(d)
        mol = mc.mol
        ncore = mc.ncore
        nact = rdm1.shape[0]
        nmo = mo_coeff.shape[0]
        Dd = numpy.dot(rdm1, d)
        Dd_full = numpy.zeros((nmo, nmo))
        Dd_full[ncore : ncore + nact, ncore : ncore + nact] = Dd

        if _opt_impl:
            Dd_full_AO = reduce(numpy.dot, (mo_coeff, Dd_full, mo_coeff.T))
            K_ao = mc._scf.get_k(dm=Dd_full_AO)
            K = reduce(numpy.dot, (mo_coeff.T, K_ao, mo_coeff))
        else:
            int2e = pyscf.ao2mo.general(
                mol,
                (mo_coeff, mo_coeff, mo_coeff, mo_coeff),
                aosym="1",
                compact=False,
            ).reshape(nmo, nmo, nmo, nmo)
            K = numpy.einsum(
                "prsq,rs->pq",
                int2e[:, ncore : ncore + nact, ncore : ncore + nact, :],
                Dd,
            )
        g3 = -0.5 * reduce(numpy.dot, (Dd_full, K, Dd_full))
        # print(g3)
        # exit(1)

        if only_g3:
            return g3
        else:
            return reduce(numpy.dot, (mo_coeff.conj().T, fock_ao, mo_coeff)) + g3


############## integral MRPT2 ##############


def fcidump_mrpt2(mol, scf, mo_coeff, nfzc, nact, nvir, filename="FCIDUMP", tol=1e-10):
    nmo = nfzc + nact + nvir
    assert nmo <= mol.nao
    # nmo = my_scf.mo_coeff.shape[1]
    nelec = mol.nelectron
    ms = 0
    # tol = 1e-10
    nuc = mol.get_enuc()
    float_format = tools.fcidump.DEFAULT_FLOAT_FORMAT

    h1e = reduce(numpy.dot, (mo_coeff.T, scf.get_hcore(), mo_coeff))
    h1e = h1e[:nmo, :nmo]

    # print(h1e)

    int2e_full = pyscf.ao2mo.full(eri_or_mol=mol, mo_coeff=mo_coeff[:, :nmo], aosym="4")
    int2e_full = pyscf.ao2mo.restore(8, int2e_full.copy(), nmo)

    OrbSym = pyscf.symm.label_orb_symm(
        mol, mol.irrep_name, mol.symm_orb, mo_coeff[:, :nmo]
    )
    OrbSymID = [pyscf.symm.irrep_name2id(mol.groupname, x) for x in OrbSym]

    with open(filename, "w") as fout:  # 8-fold symmetry
        tools.fcidump.write_head(fout, nmo, nelec, ms, OrbSymID)
        output_format = float_format + " %4d %4d %4d %4d\n"
        for p in range(nmo):

            ncore_p = 0
            nvirt_p = 0
            if p < nfzc:
                ncore_p += 1
            if p >= (nfzc + nact):
                nvirt_p += 1

            for q in range(p + 1):

                ncore_q = ncore_p
                nvirt_q = nvirt_p
                if q < nfzc:
                    ncore_q += 1
                if q >= (nfzc + nact):
                    nvirt_q += 1

                for r in range(p + 1):

                    ncore_r = ncore_q
                    nvirt_r = nvirt_q
                    if r < nfzc:
                        ncore_r += 1
                    if r >= (nfzc + nact):
                        nvirt_r += 1

                    if (p == q) or (p == r) or (q == r):
                        if p > r:
                            for s in range(r + 1):
                                indx = _combine4(p, q, r, s)
                                if abs(int2e_full[indx]) > tol:
                                    fout.write(
                                        output_format
                                        % (int2e_full[indx], p + 1, q + 1, r + 1, s + 1)
                                    )
                        else:
                            for s in range(q + 1):
                                indx = _combine4(p, q, r, s)
                                if abs(int2e_full[indx]) > tol:
                                    fout.write(
                                        output_format
                                        % (int2e_full[indx], p + 1, q + 1, r + 1, s + 1)
                                    )
                    else:

                        if (ncore_r > 2) or (nvirt_r > 2):
                            continue

                        if p > r:
                            for s in range(r + 1):
                                ncore_s = ncore_r
                                nvirt_s = nvirt_r
                                if s < nfzc:
                                    ncore_s += 1
                                if s >= (nfzc + nact):
                                    nvirt_s += 1
                                indx = _combine4(p, q, r, s)
                                if (p == s) or (q == s) or (r == s):
                                    if abs(int2e_full[indx]) > tol:
                                        fout.write(
                                            output_format
                                            % (
                                                int2e_full[indx],
                                                p + 1,
                                                q + 1,
                                                r + 1,
                                                s + 1,
                                            )
                                        )
                                else:
                                    if (ncore_s <= 2) and (nvirt_s <= 2):
                                        if abs(int2e_full[indx]) > tol:
                                            fout.write(
                                                output_format
                                                % (
                                                    int2e_full[indx],
                                                    p + 1,
                                                    q + 1,
                                                    r + 1,
                                                    s + 1,
                                                )
                                            )
                                    else:
                                        if nvirt_s > 2:
                                            break
                        else:
                            for s in range(q + 1):
                                ncore_s = ncore_r
                                nvirt_s = nvirt_r
                                if s < nfzc:
                                    ncore_s += 1
                                if s >= (nfzc + nact):
                                    nvirt_s += 1
                                indx = _combine4(p, q, r, s)
                                if (p == s) or (q == s) or (r == s):
                                    if abs(int2e_full[indx]) > tol:
                                        fout.write(
                                            output_format
                                            % (
                                                int2e_full[indx],
                                                p + 1,
                                                q + 1,
                                                r + 1,
                                                s + 1,
                                            )
                                        )
                                else:
                                    if (ncore_s <= 2) and (nvirt_s <= 2):
                                        if abs(int2e_full[indx]) > tol:
                                            fout.write(
                                                output_format
                                                % (
                                                    int2e_full[indx],
                                                    p + 1,
                                                    q + 1,
                                                    r + 1,
                                                    s + 1,
                                                )
                                            )
                                    else:
                                        if nvirt_s > 2:
                                            break

        tools.fcidump.write_hcore(fout, h1e, nmo, tol=tol, float_format=float_format)
        output_format = float_format + "  0  0  0  0\n"
        fout.write(output_format % nuc)


############## integral MRPT2 outcore ##############


def fcidump_mrpt2_outcore(
    mol, scf, mo_coeff, nfzc, nact, nvir, filename="FCIDUMP", tol=1e-8
):
    nmo = nfzc + nact + nvir
    assert nmo <= mol.nao
    # nmo = my_scf.mo_coeff.shape[1]
    nelec = mol.nelectron
    ms = 0
    # tol = 1e-10
    nuc = mol.get_enuc()
    float_format = tools.fcidump.DEFAULT_FLOAT_FORMAT

    h1e = reduce(numpy.dot, (mo_coeff.T, scf.get_hcore(), mo_coeff))
    h1e = h1e[:nmo, :nmo]

    ftmp = tempfile.NamedTemporaryFile()
    print('MO integrals are saved in file  %s  under dataset "eri_mo"' % ftmp.name)
    pyscf.ao2mo.kernel(mol, mo_coeff[:, :nmo], ftmp.name)

    OrbSym = pyscf.symm.label_orb_symm(
        mol, mol.irrep_name, mol.symm_orb, mo_coeff[:, :nmo]
    )
    OrbSymID = [pyscf.symm.irrep_name2id(mol.groupname, x) for x in OrbSym]

    # dump

    with open(filename, "w") as fout:  # 8-fold symmetry
        tools.fcidump.write_head(fout, nmo, nelec, ms, OrbSymID)
        output_format = float_format + " %4d %4d %4d %4d\n"

        # eri

        with h5py.File(ftmp.name) as eri_file:
            for p in range(nmo):

                ncore_p = 0
                nvirt_p = 0
                if p < nfzc:
                    ncore_p += 1
                if p >= (nfzc + nact):
                    nvirt_p += 1

                for q in range(p + 1):

                    ncore_q = ncore_p
                    nvirt_q = nvirt_p
                    if q < nfzc:
                        ncore_q += 1
                    if q >= (nfzc + nact):
                        nvirt_q += 1

                    eri_pq = eri_file["eri_mo"][_combine2(p, q)]

                    # print("shape %d %d " % (p, q), eri_pq.shape)

                    for r in range(p + 1):

                        ncore_r = ncore_q
                        nvirt_r = nvirt_q
                        if r < nfzc:
                            ncore_r += 1
                        if r >= (nfzc + nact):
                            nvirt_r += 1

                        if (p == q) or (p == r) or (q == r):
                            if p > r:
                                for s in range(r + 1):
                                    indx = _combine2(r, s)
                                    if abs(eri_pq[indx]) > tol:
                                        fout.write(
                                            output_format
                                            % (eri_pq[indx], p + 1, q + 1, r + 1, s + 1)
                                        )
                            else:
                                for s in range(q + 1):
                                    indx = _combine2(r, s)
                                    if abs(eri_pq[indx]) > tol:
                                        fout.write(
                                            output_format
                                            % (eri_pq[indx], p + 1, q + 1, r + 1, s + 1)
                                        )
                        else:

                            if (ncore_r > 2) or (nvirt_r > 2):
                                continue

                            if p > r:
                                for s in range(r + 1):
                                    ncore_s = ncore_r
                                    nvirt_s = nvirt_r
                                    if s < nfzc:
                                        ncore_s += 1
                                    if s >= (nfzc + nact):
                                        nvirt_s += 1
                                    indx = _combine2(r, s)
                                    if (p == s) or (q == s) or (r == s):
                                        if abs(eri_pq[indx]) > tol:
                                            fout.write(
                                                output_format
                                                % (
                                                    eri_pq[indx],
                                                    p + 1,
                                                    q + 1,
                                                    r + 1,
                                                    s + 1,
                                                )
                                            )
                                    else:
                                        if (ncore_s <= 2) and (nvirt_s <= 2):
                                            if abs(eri_pq[indx]) > tol:
                                                fout.write(
                                                    output_format
                                                    % (
                                                        eri_pq[indx],
                                                        p + 1,
                                                        q + 1,
                                                        r + 1,
                                                        s + 1,
                                                    )
                                                )
                                        else:
                                            if nvirt_s > 2:
                                                break
                            else:
                                for s in range(q + 1):
                                    ncore_s = ncore_r
                                    nvirt_s = nvirt_r
                                    if s < nfzc:
                                        ncore_s += 1
                                    if s >= (nfzc + nact):
                                        nvirt_s += 1
                                    indx = _combine2(r, s)
                                    if (p == s) or (q == s) or (r == s):
                                        if abs(eri_pq[indx]) > tol:
                                            fout.write(
                                                output_format
                                                % (
                                                    eri_pq[indx],
                                                    p + 1,
                                                    q + 1,
                                                    r + 1,
                                                    s + 1,
                                                )
                                            )
                                    else:
                                        if (ncore_s <= 2) and (nvirt_s <= 2):
                                            if abs(eri_pq[indx]) > tol:
                                                fout.write(
                                                    output_format
                                                    % (
                                                        eri_pq[indx],
                                                        p + 1,
                                                        q + 1,
                                                        r + 1,
                                                        s + 1,
                                                    )
                                                )
                                        else:
                                            if nvirt_s > 2:
                                                break

        # h1e and nuc

        tools.fcidump.write_hcore(fout, h1e, nmo, tol=tol, float_format=float_format)
        output_format = float_format + "  0  0  0  0\n"
        fout.write(output_format % nuc)


if __name__ == "__main__":

    from pyscf_util.Integrals.integral_CASCI import dump_heff_casci
    from pyscf_util.MeanField import iciscf
    from pyscf_util.iCIPT2.iCIPT2 import kernel
    import os
    from pyscf.tools import fcidump
    from pyscf import symm
    from pyscf_util.Integrals.integral_sfX2C import fcidump_sfx2c
    from pyscf_util.File.file_cmoao import Dump_Cmoao

    def OrbSymInfo(Mol, mo_coeff):
        IRREP_MAP = {}
        nsym = len(Mol.irrep_name)
        for i in range(nsym):
            IRREP_MAP[Mol.irrep_name[i]] = i
        # print(IRREP_MAP)

        OrbSym = pyscf.symm.label_orb_symm(Mol, Mol.irrep_name, Mol.symm_orb, mo_coeff)
        IrrepOrb = []
        for i in range(len(OrbSym)):
            IrrepOrb.append(symm.irrep_name2id(Mol.groupname, OrbSym[i]))
        return IrrepOrb

    ### take Cr2 as the show case ###

    Mol = pyscf.gto.Mole()
    Mol.atom = """
Cr     0.0000      0.0000   %f 
Cr     0.0000      0.0000  -%f 
""" % (
        1.68 / 2,
        1.68 / 2,
    )
    Mol.basis = "def2-svp"
    Mol.symmetry = "Dooh"
    Mol.spin = 2
    Mol.charge = 0
    Mol.verbose = 4
    Mol.unit = "angstorm"
    Mol.build()
    SCF = pyscf.scf.sfx2c(pyscf.scf.RHF(Mol))
    SCF.max_cycle = 32
    SCF.conv_tol = 1e-9
    SCF.run()

    Mol.spin = 0
    Mol.build()

    norb = 12
    nelec = 12
    CASSCF_Driver = pyscf.mcscf.CASSCF(SCF, norb, nelec)

    cas_space_symmetry = {
        "A1u": 2,  # 5
        "A1g": 2,  # 0
        "E1ux": 1,  # 7
        "E1gy": 1,  # 3
        "E1gx": 1,  # 2
        "E1uy": 1,  # 6
        "E2gy": 1,  # 1
        "E2gx": 1,  # 0
        "E2uy": 1,  # 4
        "E2ux": 1,  # 5
    }
    ### generate init guess for CAS ###

    mo_init = pyscf.mcscf.sort_mo_by_irrep(
        CASSCF_Driver, CASSCF_Driver.mo_coeff, cas_space_symmetry
    )  # right!
    SCF.mo_coeff = mo_init
    CASSCF_Driver = iciscf.iCISCF(SCF, norb, nelec, cmin=0.0)
    CASSCF_Driver.canonicalization = True

    ### run ###

    energy, _, _, mo_coeff, mo_energy = CASSCF_Driver.kernel(mo_coeff=mo_init)

    Mol.symmetry = "D2h"
    Mol.build()

    ### call icipt2 to generate rdm1 ###

    # mo_coeff = CASSCF_Driver.mo_coeff

    dump_heff_casci(
        Mol,
        CASSCF_Driver,
        mo_coeff[:, :18],
        mo_coeff[:, 18:30],
        _filename="FCIDUMP_Cr2",
    )

    kernel(
        IsCSF=True,
        task_name="cr2_rdm1",
        fcidump="FCIDUMP_Cr2",
        segment="0 0 6 6 0 0",
        nelec_val=12,
        rotatemo=0,
        cmin=0.0,
        perturbation=0,
        dumprdm=1,
        relative=0,
        Task="0 0 1 1",
        inputocfg=0,
        etol=1e-10,
        selection=1,
        doublegroup=None,
        direct=None,
        start_with=None,
        end_with=[".csv"],
    )

    os.system("mv rdm1.csv cr2_rdm1.csv")

    mo_coeff = CASSCF_Driver.mo_coeff
    rdm1 = file_rdm.ReadIn_rdm1("cr2_rdm1", 12, 12)

    print(rdm1)
    print(get_generalized_fock(CASSCF_Driver, mo_coeff, rdm1))
    print(mo_energy)

    gfock = get_generalized_fock(CASSCF_Driver, mo_coeff, rdm1)

    g3_1 = get_generalized_fock(CASSCF_Driver, mo_coeff, rdm1, True, True, True)
    g3_2 = get_generalized_fock(CASSCF_Driver, mo_coeff, rdm1, True, True, False)
    assert numpy.allclose(g3_1, g3_2)

    # fcidump #

    orbsym = OrbSymInfo(Mol, CASSCF_Driver.mo_coeff)
    # fcidump.from_mo(
    #     Mol, "FCIDUMP_Cr2_Benchmark", CASSCF_Driver.mo_coeff, orbsym, tol=1e-12
    # )
    fcidump_sfx2c(Mol, SCF, mo_coeff, "FCIDUMP_Cr2_Benchmark", 1e-10)
    fcidump_mrpt2_outcore(
        Mol, SCF, mo_coeff, 18, 12, Mol.nao - 30, "FCIDUMP_Cr2_outcore", 1e-10
    )
    fcidump_mrpt2(Mol, SCF, mo_coeff, 18, 12, Mol.nao - 30, "FCIDUMP_Cr2_incore", 1e-10)

    Dump_Cmoao("gfock", gfock)
