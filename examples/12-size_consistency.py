from pyscf import gto, scf, symm

# from pyscf_util.iCIPT2.iCIPT2_coov import kernel
from pyscf import tools

# from pyscf_util.iCIPT2.iCIPT2 import kernel
import pyscf
from pyscf_util.Integrals.integral_MRPT2 import get_generalized_fock
from pyscf_util.MeanField.iciscf import iCI
from pyscf_util.Integrals.integral_CASCI import dump_heff_casci
from pyscf_util.iCIPT2.iCIPT2 import kernel

import numpy as np


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


mol = gto.M(
    verbose=4,
    atom="""
            C   0.000000000000       0.000000000000      -0.621265
            C   0.000000000000       0.000000000000       0.621265
            """,
    basis={"C": "cc-pvdz", "O": "cc-pvdz"},
    spin=0,
    charge=4,
    symmetry="d2h",
)
mol.build()
mf = scf.RHF(mol)
mf.kernel()

mf.analyze()

# exit(1)

norb = 8
nelec = 4
CASSCF_Driver = pyscf.mcscf.CASSCF(mf, norb, nelec)
CASSCF_Driver.fcisolver = iCI(
    mol=mol,
    cmin=0.0,
    state=[[0, 0, 1]],
    tol=1e-12,
    mo_coeff=mf.mo_coeff,
    taskname="iCI0",
)
CASSCF_Driver.mc1step()

### dump heff and generate gfock ###

mo_coeff = CASSCF_Driver.mo_coeff

# start, end = 10, 28  # 包含第2列到第6列
# cols_to_shuffle = mo_coeff[:, start:end]
# n_cols = end - start
# perm = np.random.permutation(n_cols)
# shuffled_cols = cols_to_shuffle[:, perm]
# mo_coeff[:, start:end] = shuffled_cols

dump_heff_casci(
    mol,
    CASSCF_Driver,
    mo_coeff[:, :2],
    mo_coeff[:, 2:10],
    _filename="FCIDUMP_C2",
)

kernel(
    IsCSF=True,
    task_name="c2_rdm1",
    fcidump="FCIDUMP_C2",
    segment="0 0 4 4 0 0",
    nelec_val=4,
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

import os
from pyscf_util.File import file_rdm, file_cmoao

os.system("mv rdm1.csv c2_rdm1.csv")

# mo_coeff = CASSCF_Driver.mo_coeff

rdm1 = file_rdm.ReadIn_rdm1("c2_rdm1", 8, 8)

gfock = get_generalized_fock(CASSCF_Driver, mo_coeff, rdm1)

from pyscf_util.mrpt2.canonicalize import *

mo_coeff_new, new_fock, orb_ene = Canonicalize(mol, mo_coeff, gfock, 2, 8, mol.nao - 10)

NVIR = 18
mo_coeff = mo_coeff_new
gfock = new_fock

gfock = gfock[: 10 + NVIR, : 10 + NVIR]

file_cmoao.Dump_Cmoao("gfock", gfock)
mf.mo_coeff = mo_coeff


mo_coeff = mo_coeff[:, : 2 + 8 + NVIR]

orb_id = OrbSymInfo(mol, mo_coeff)

h1e_mol, h2e_mol, _, _, _, _ = dump_heff_casci(
    mol,
    CASSCF_Driver,
    mo_coeff[:, :0],
    mo_coeff[:, :],
    None,
)

dump_heff_casci(
    mol,
    CASSCF_Driver,
    mo_coeff[:, :0],
    mo_coeff,
    "FCIDUMP_C2",
)

from pyscf_util.mrpt2.un_nevpt2 import *

h1e_dyall, h2e_dyall, ecore_dyall, _, _, _ = fcidump_Dyall(
    mol,
    CASSCF_Driver,
    mo_coeff,
    np.diag(gfock),
    2,
    8,
    NVIR,
    _filename=None,
)

fcidump_Dyall(
    mol,
    CASSCF_Driver,
    mo_coeff,
    np.diag(gfock),
    2,
    8,
    NVIR,
    _filename="FCIDUMP_C2_Dyall",
)

from pyscf_util.mrpt2.size_consistency import *

gfock1 = gfock_size_consistency(
    # mol 1 #
    gfock,
    2,
    8,
    NVIR,
    # mol 2 #
    gfock,
    2,
    8,
    NVIR,
)

print(np.diag(gfock[:2, :2]))
print(np.diag(gfock1[:4, :4]))
# print(np.diag(gfock[10:,10:]))
# print(np.diag(gfock1[20:,20:]))

file_cmoao.Dump_Cmoao("gfock2", gfock1)

ecore_mol = mf.energy_nuc()

fcidump_size_consistency(
    # mol 1 #
    h1e_mol,
    h2e_mol,
    ecore_mol,
    orb_id,
    2,
    8,
    NVIR,
    8,
    # mol2 #
    h1e_mol,
    h2e_mol,
    ecore_mol,
    orb_id,
    2,
    8,
    NVIR,
    8,
    # filename #
    f"FCIDUMP_C2_C2",
)

fcidump_size_consistency(
    # mol 1 #
    h1e_dyall,
    h2e_dyall,
    ecore_dyall,
    orb_id,
    2,
    8,
    NVIR,
    8,
    # mol2 #
    h1e_dyall,
    h2e_dyall,
    ecore_dyall,
    orb_id,
    2,
    8,
    NVIR,
    8,
    # filename #
    f"FCIDUMP_C2_C2_Dyall",
)
