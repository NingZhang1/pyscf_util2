from pyscf import gto, scf
from pyscf import tools

from pyscf_util._atmMinCASOrb._orb_loader import LoadAtmHFOrb
from pyscf_util.misc.mole import get_mol
from pyscf_util.MeanField.scf import kernel as scf_kernel
from pyscf_util.Relativisitc.sfX2C_soDKH import fetch_X2C_soDKH1
import pyscf_util.File.file_sodkh13 as file_sodkh13
from pyscf_util.iCIPT2.iCIPT2 import kernel
from pyscf import data
import numpy

Cl_atm = get_mol("Cl 0 0 0", 0, 1, basis="unc-cc-pvdz-dk", symmetry="d2h", verbose=11)
Cl_scf = scf_kernel(Cl_atm, sfx1e=True, run=True)
Cl_res = LoadAtmHFOrb("Cl", 0, "unc-cc-pvdz-dk", with_sfx2c=True, rerun=True)

Cl_scf.mo_coeff = Cl_res["mo_coeff"]
Cl_scf.mo_energy = Cl_res["mo_energy"]

# dump FCIDUMP

FCIDUMP_NAME = "FCIDUMP"
RELDUMP_NAME = "RELDUMP"

tools.fcidump.from_scf(Cl_scf, FCIDUMP_NAME, 1e-10)

# dm1 = Cl_scf.make_rdm1()  # TODO: check X2C !!
dm1 = numpy.zeros((Cl_atm.nao, Cl_atm.nao))
nelectrons = Cl_atm.nelectron
print("nelectrons = ", nelectrons)
print("core = ", (nelectrons - 5) // 2)
for i in range((nelectrons - 5) // 2):
    dm1[i, i] = 2
for i in range((nelectrons - 5) // 2, (nelectrons - 5) // 2 + 3):
    dm1[i, i] = 5.0 / 3.0
# print(dm1[:9, :9])

# fetch X2C integrals #

hsf = numpy.zeros((1, Cl_atm.nao, Cl_atm.nao))
hso = fetch_X2C_soDKH1(
    Cl_atm, Cl_scf, Cl_res["mo_coeff"], dm1, _test=False, _get_1e=True, _get_2e=True
)
hso *= (data.nist.ALPHA**2) / 4.0
# print(hso[0][:9, :9])
hso = numpy.vstack((hso, hsf))
# print(hso.shape)
file_sodkh13.Dump_Relint_iCI("RELDUMP", hso, Cl_atm.nao)
# exit(1)

# kernel(
#     True,
#     task_name="iCIPT2_Cl",
#     fcidump=FCIDUMP_NAME,
#     segment="5 0 4 4 %d 0" % (Cl_atm.nao - 13),
#     nelec_val=7,
#     relative=1,
#     cmin="1e-4",
#     perturbation=1,
#     Task="1 5 1 1 1 6 1 1 1 7 1 1",
#     doublegroup="d2h",
# )
