from pyscf import gto, scf

# from pyscf_util.iCIPT2.iCIPT2_4C_d2h import kernel
from pyscf_util.Relativisitc.integral_4C import (
    FCIDUMP_Rela4C_SU2,
)
from pyscf_util.Relativisitc.integral_4C_opt import (
    FCIDUMP_Rela4C_SU2_opt,
)
from pyscf_util.Relativisitc.integral_4C_binary import (
    FCIDUMP_Rela4C_SU2_binary,
)

mol = gto.M(
    atom="Cl 0 0 0",
    basis="unc-cc-pvtz-dk",
    verbose=5,
    charge=-1,
    spin=0,
    symmetry="d2h",
)
mol.build()
mf = scf.dhf.RDHF(mol)
mf.conv_tol = 1e-12
mf.kernel()

mf.with_breit = True
mf.kernel()

FCIDUMP_NAME = "FCIDUMP_Cl"

FCIDUMP_Rela4C_SU2_opt(mol, mf, True, filename=FCIDUMP_NAME + "_OPT", mode="outcore")
FCIDUMP_Rela4C_SU2_binary(
    mol,
    mf,
    True,
    filename=FCIDUMP_NAME + "_head",
    binary_filename=FCIDUMP_NAME + ".bin",
    mode="outcore",
)
