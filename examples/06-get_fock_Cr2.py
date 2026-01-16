from pyscf import gto, scf

# from pyscf_util.iCIPT2.iCIPT2_coov import kernel
from pyscf import tools

# from pyscf_util.iCIPT2.iCIPT2 import kernel
import pyscf
from pyscf_util.Integrals.integral_MRPT2 import get_generalized_fock
from pyscf_util.MeanField.iciscf import iCI
from pyscf_util.Integrals.integral_CASCI import dump_heff_casci
from pyscf_util.iCIPT2.iCIPT2 import kernel
from pyscf_util.MeanField import iciscf

cas_space_symmetry = {
    "A1u": 1 + 1,  # 5
    "A1g": 1 + 1,  # 0
    "E1ux": 1,  # 7
    "E1gy": 1,  # 3
    "E1gx": 1,  # 2
    "E1uy": 1,  # 6
    "E2gy": 1,  # 1
    "E2gx": 1,  # 0
    "E2uy": 1,  # 4
    "E2ux": 1,  # 5
}

cas_space_symmetry = {
    "Ag": 3,
    "B1g": 1,
    "B2g": 1,
    "B3g": 1,
    "au": 1,
    "b1u": 3,
    "b2u": 1,
    "b3u": 1,
}


mol = gto.M(
    verbose=4,
    atom="""
            Cr   0.000000000000       0.000000000000      -0.84
            Cr   0.000000000000       0.000000000000       0.84
            """,
    basis={"Cr": "vdz", "O": "cc-pvdz"},
    spin=0,
    charge=0,
    symmetry="d2h",
    unit="angstorm",
)
mol.build()
mf = scf.RHF(mol)
mf = pyscf.scf.sfx2c(pyscf.scf.RHF(mol))
mf.kernel()

norb = 12
nelec = 12
CASSCF_Driver = pyscf.mcscf.CASSCF(mf, norb, nelec)
mo_init = pyscf.mcscf.sort_mo_by_irrep(
    CASSCF_Driver, CASSCF_Driver.mo_coeff, cas_space_symmetry
)  # right!
mf.mo_coeff = mo_init
# CASSCF_Driver = pyscf.mcscf.CASSCF(mf, norb, nelec)
# solver1 = pyscf.fci.direct_spin1_symm.FCI(mol)
# solver1.wfnsym = "ag"
# solver1.nroots = 1
# solver1.spin = 0
# CASSCF_Driver.mc1step()
# CASSCF_Driver.fcisolver = iCI(
#     mol=mol,
#     cmin=0.0,
#     state=[[0, 0, 1]],
#     tol=1e-12,
#     mo_coeff=mf.mo_coeff,
#     taskname="iCI0",
# )
# CASSCF_Driver.mc1step()
CASSCF_Driver = iciscf.iCISCF(mf, norb, nelec, cmin=0.0)
energy, _, _, mo_coeff, _ = CASSCF_Driver.kernel(mo_coeff=mo_init)

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
file_cmoao.Dump_Cmoao("gfock", gfock)
gfock = get_generalized_fock(CASSCF_Driver, mo_coeff, rdm1, True)
file_cmoao.Dump_Cmoao("gfock2", gfock)
mf.mo_coeff = mo_coeff

tools.fcidump.from_scf(mf, "FCIDUMP_Cr2", 1e-10)


# MRPT2 2

mc = pyscf.mcscf.CASSCF(mf, norb, nelec)
solver1 = pyscf.fci.direct_spin1_symm.FCI(mol)
solver1.wfnsym = "ag"
solver1.nroots = 1
solver1.spin = 0
mc.mc1step()

from pyscf import mrpt

# mrpt.nevpt2.sc_nevpt(mc)
