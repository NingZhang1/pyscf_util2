import pyscf
from pyscf import gto, scf


def get_orbsym(mol, mocoeff):

    OrbSym = pyscf.symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, mocoeff)
    OrbSymID = [pyscf.symm.irrep_name2id(mol.groupname, x) for x in OrbSym]

    return OrbSymID, OrbSym


def get_mol(
    xyz, charge=0, spin=0, basis="6-31G(d)", symmetry="", verbose=4, unit="angstorm"
):
    mol = pyscf.gto.M(
        verbose=verbose,
        atom=xyz,
        basis=basis,
        spin=spin,
        charge=charge,
        symmetry=symmetry,
        unit=unit,
    )
    mol.build()
    return mol


###################
# fake mol #
###################


def get_fake_mol(nelectron, nao, verbose=4, ao_labels=None):
    mol_fake = gto.M(verbose=verbose, dump_input=False)
    mol_fake.nelectron = int(nelectron)
    mol_fake.nao_nr = lambda *args: nao
    mol_fake.incore_anyway = True
    if ao_labels is not None:
        mol_fake.ao_labels = lambda *args: ao_labels

    return mol_fake


###################
# fake mf #
###################

import numpy as np


def get_fake_mf(mol, h0, hcore, eri, ovlp=None, run=True):
    if ovlp is None:
        ovlp = np.identity(mol.nao_nr())
    mf = scf.RHF(mol)
    mf.energy_nuc = lambda *args: h0
    mf.get_hcore = lambda *args: hcore
    mf.get_ovlp = lambda *args: ovlp
    mf._eri = eri
    if run:
        mf.kernel()
    return mf


if __name__ == "__main__":

    Mol = pyscf.gto.Mole()
    Mol.atom = """
    C      0.0000      1.396792    0.0000
    C      0.0000     -1.396792    0.0000
    C      1.209657    0.698396    0.0000
    C     -1.209657   -0.698396    0.0000
    C     -1.209657    0.698396    0.0000
    C      1.209657   -0.698396    0.0000
    H      0.0000      2.484212    0.0000
    H      2.151390    1.242106    0.0000
    H     -2.151390   -1.242106    0.0000
    H     -2.151390    1.242106    0.0000
    H      2.151390   -1.242106    0.0000
    H      0.0000     -2.484212    0.0000
    """
    Mol.basis = "ccpvdz"
    Mol.symmetry = "D2h"
    Mol.spin = 0
    Mol.verbose = 4
    Mol.build()

    SCF = pyscf.scf.RHF(Mol)
    SCF.run()

    SCF.analyze()

    nact_mol = 6
    nelec_act = 6

    cas_space_symmetry = {
        "Au": 1,  # 5
        "B1u": 2,  # 0
        "B2g": 1,  # 7
        "B3g": 2,  # 3
    }
    core_space_symmetry = {
        "Ag": 6,  # 5
        "B1g": 3,  # 0
        "B2g": 0,  # 7
        "B3g": 0,  # 3
        "Au": 0,  # 5
        "B1u": 0,  # 0
        "B2u": 5,  # 7
        "B3u": 4,  # 3
    }

    from pyscf import mcscf
    from pyscf import mrpt
    from pyscf_util.iCIPT2.iCIPT2 import kernel

    CASSCF_Driver = pyscf.mcscf.CASSCF(SCF, nact_mol, nelec_act)
    mo_init = pyscf.mcscf.sort_mo_by_irrep(
        CASSCF_Driver, CASSCF_Driver.mo_coeff, cas_space_symmetry, core_space_symmetry
    )  # right!
    SCF.mo_coeff = mo_init.copy()

    CASSCF_Driver = pyscf.mcscf.CASSCF(SCF, nact_mol, nelec_act)
    energy, _, _, mo_coeff, _ = CASSCF_Driver.kernel(mo_coeff=mo_init)

    mrpt.nevpt2.sc_nevpt(CASSCF_Driver)

    # fake mol and fake mf #

    from pyscf_util.Integrals.integral_CASCI import *

    ecore_mol = Mol.energy_nuc()
    h1e_mol, h2e_mol, _, _, _, _ = dump_heff_casci(
        Mol,
        CASSCF_Driver,
        mo_init[:, :0],
        mo_init,
        None,
    )

    fake_mol = get_fake_mol(Mol.nelectron, Mol.nao, 4, Mol.ao_labels())
    fake_mf = get_fake_mf(fake_mol, ecore_mol, h1e_mol, h2e_mol, run=False)
    fake_mf.mo_coeff = np.identity(fake_mol.nao_nr())
    fake_casscf = pyscf.mcscf.CASSCF(fake_mf, nact_mol, nelec_act)
    energy, _, _, mo_coeff, _ = fake_casscf.kernel(
        mo_coeff=np.identity(fake_mol.nao_nr())
    )
    # fake_casscf.mc1step()

    mrpt.nevpt2.sc_nevpt(fake_casscf)

    # pc-nevpt2 #

    mo_coeff = fake_casscf.mo_coeff
    ncore = fake_casscf.ncore

    dump_heff_casci(
        fake_mol,
        fake_casscf,
        mo_coeff[:, :ncore],
        mo_coeff[:, ncore : ncore + nact_mol],
        _filename="FCIDUMP_benzene",
    )

    kernel(
        IsCSF=True,
        task_name="benzene_rdm12",
        fcidump="FCIDUMP_benzene",
        segment="0 0 3 3 0 0",
        nelec_val=6,
        rotatemo=0,
        cmin=0.0,
        perturbation=0,
        dumprdm=2,
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

    # run iCI to dump rdm1 and rdm2

    file_rdm1 = "rdm1.csv"
    file_rdm2 = "rdm2.csv"

    from pyscf_util.mrpt2.nevpt2 import *

    # run iCIPT2-NEVPT2

    res, eris = sc_nevpt2_ici(fake_casscf, None, file_rdm1, file_rdm2)

    # print res

    for k, v in res.items():
        print(k)
        print("norm %15.8f e %15.8f" % (v["norm"], v["e"]))
        print("-" * 100)

    res, eris = pc_nevpt2_ici(fake_casscf, eris, file_rdm1, file_rdm2)

    for k, v in res.items():
        print(k)
        print("e %15.8f" % (v["e"]))
        print("-" * 100)
