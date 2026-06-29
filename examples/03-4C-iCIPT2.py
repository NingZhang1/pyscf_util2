from pyscf import gto, scf

# from pyscf_util.iCIPT2.iCIPT2_4C_d2h import kernel
from pyscf_util.Relativisitc.integral_4C import (
    FCIDUMP_Rela4C_SU2,
)

mol = gto.M(
    atom="F 0 0 0",
    basis="unc-cc-pvdz-dk",
    verbose=5,
    charge=-1,
    spin=0,
    symmetry="d2h",
)
mol.build()
mf = scf.dhf.RDHF(mol)
mf.conv_tol = 1e-12
mf.kernel()

FCIDUMP_NAME = "FCIDUMP_F_Coulomb"
FCIDUMP_Rela4C_SU2(
    mol, mf, False, filename=FCIDUMP_NAME, mode="outcore", debug=True, npes=18
)
# FCIDUMP_NAME = "FCIDUMP_F_Coulomb_no2e"
# FCIDUMP_Rela4C_SU2(
#     mol, mf, False, filename=FCIDUMP_NAME, mode="outcore", debug=True, no_2e=True
# )

mf.with_breit = True
mf.kernel()

FCIDUMP_NAME = "FCIDUMP_F"
FCIDUMP_Rela4C_SU2(
    mol, mf, True, filename=FCIDUMP_NAME, mode="outcore", debug=True, npes=18
)

# FCIDUMP_NAME = "FCIDUMP_F_no2e"
# FCIDUMP_Rela4C_SU2(
#     mol, mf, True, filename=FCIDUMP_NAME, mode="outcore", debug=True, no_2e=True
# )

# kernel(
#     True,
#     task_name="iCIPT2_4C_F",
#     fcidump=FCIDUMP_NAME,
#     segment="1 0 4 4 %d 0" % (mol.nao - 9),
#     nelec_val=7,
#     cmin="1e-6",
#     perturbation=1,
#     Task="1 1 3 1 1 1",
#     end_with=".PrimeSpace",
# )

mol = gto.M(
    atom="O 0 0 0",
    basis="unc-cc-pvdz-dk",
    verbose=5,
    charge=-2,
    spin=0,
    symmetry="d2h",
)
mol.build()
mf = scf.dhf.RDHF(mol)
mf.conv_tol = 1e-12
mf.kernel()

FCIDUMP_NAME = "FCIDUMP_O_Coulomb"
FCIDUMP_Rela4C_SU2(
    mol, mf, False, filename=FCIDUMP_NAME, mode="outcore", debug=True, npes=18
)

mf.with_breit = True
mf.kernel()

FCIDUMP_NAME = "FCIDUMP_O"
FCIDUMP_Rela4C_SU2(
    mol, mf, True, filename=FCIDUMP_NAME, mode="outcore", debug=True, npes=18
)
