from pyscf import gto, scf

# from pyscf_util.iCIPT2.iCIPT2_4C_d2h import kernel
from pyscf_util.Relativisitc.integral_4C import (
    FCIDUMP_Rela4C,
)
from pyscf_util.Relativisitc.double_group import time_reversal_symm_adapted

mol = gto.M(
    atom="""
    Be 0 0 0
    Be 3 0 0 
    """,
    basis="sto-3g",
    verbose=5,
    charge=0,
    spin=0,
    symmetry="d2h",
)
mol.build()
mf = scf.dhf.RDHF(mol)
mf.conv_tol = 1e-12
mf.kernel()

# mo_coeff = time_reversal_symm_adapted(mol, mf.mo_coeff)

FCIDUMP_NAME = "FCIDUMP_Be2_Coulomb"
FCIDUMP_Rela4C(mol, mf, False, filename=FCIDUMP_NAME, mode="outcore", debug=True)

# FCIDUMP_NAME = "FCIDUMP_Be2_Coulomb_no2e"
# FCIDUMP_Rela4C(mol, mf, False, filename=FCIDUMP_NAME, mode="outcore", no_2e=True, debug=True)

mf.with_breit = True
mf.kernel()

FCIDUMP_NAME = "FCIDUMP_Be2"
FCIDUMP_Rela4C(mol, mf, True, filename=FCIDUMP_NAME, mode="outcore", debug=True)
