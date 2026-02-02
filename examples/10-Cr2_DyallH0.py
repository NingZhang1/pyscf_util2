from pyscf import gto, scf

# from pyscf_util.iCIPT2.iCIPT2_coov import kernel
from pyscf import tools

# from pyscf_util.iCIPT2.iCIPT2 import kernel
import pyscf
from pyscf_util.Integrals.integral_MRPT2 import get_generalized_fock
from pyscf_util.MeanField.iciscf import iCI
from pyscf_util.Integrals.integral_CASCI import dump_heff_casci
from pyscf_util.iCIPT2.iCIPT2 import kernel
from pyscf_util.mrpt2.canonicalize import *
from pyscf_util.mrpt2.un_nevpt2 import *

mol = gto.M(
    verbose=10,
    atom="""
            Cr   0.000000000000       0.000000000000      -0.84
            Cr   0.000000000000       0.000000000000       0.84
            """,
    basis={"Cr": "cc-pvdz-dk", "O": "cc-pvdz-dk"},
    spin=12,
    charge=0,
    symmetry="d2h",
)
mol.build()
mf = pyscf.scf.sfx2c(scf.RHF(mol))
mf.kernel()

norb = 12
nelec = 12
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
CASSCF_Driver.analyze()

### dump heff and generate gfock ###

mo_coeff = CASSCF_Driver.mo_coeff

dump_heff_casci(
    mol,
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

import os
from pyscf_util.File import file_rdm, file_cmoao

os.system("mv rdm1.csv cr2_rdm1.csv")

mo_coeff = CASSCF_Driver.mo_coeff
rdm1 = file_rdm.ReadIn_rdm1("cr2_rdm1", 12, 12)

gfock = get_generalized_fock(CASSCF_Driver, mo_coeff, rdm1)

### 正则化 ###

mo_coeff_new, new_fock, orb_ene = Canonicalize(
    mol, mo_coeff, gfock, 18, 12, mol.nao - 30
)

for i in range(mol.nao):
    for j in range(mol.nao):
        if abs(new_fock[i, j]) > 1e-8:
            print("Fock element:", i, j, new_fock[i, j])

for i in range(mol.nao):
    print("Orbital energy:", i, orb_ene[i])

file_cmoao.Dump_Cmoao("Cr2_gfock", new_fock[10:42, 10:42])

# exit(1)

### dump fcidump ###

dump_heff_casci(
    mol,
    CASSCF_Driver,
    mo_coeff_new[:, :10],
    mo_coeff_new[:, 10 : 18 + 12 + 4],
    _filename="FCIDUMP_Cr2",
)

fcidump_Dyall(
    mol,
    CASSCF_Driver,
    mo_coeff_new,
    orb_ene,
    8,
    12,
    # mol.nao - 30,
    4,
    10,
    _filename="FCIDUMP_Dyall_Cr2",
)
