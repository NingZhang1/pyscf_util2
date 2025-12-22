from pyscf import gto, scf

from pyscf import tools

# from pyscf_util.iCIPT2.iCIPT2 import kernel
import pyscf

mol = gto.M(
    verbose=10,
    atom="""
            Ti   0.000000000000       0.000000000000      0.000000000000
            """,
    basis={"Ti": "cc-pvdz", "O": "cc-pvdz"},
    spin=0,
    charge=4,
    symmetry="d2h",
)
mol.build()
mf = scf.RHF(mol)
mf.kernel()
mf.analyze()

FCIDUMP_NAME = "FCIDUMP_Ti_Atm"

tools.fcidump.from_scf(mf, FCIDUMP_NAME, 1e-10)

# macro 1 1 3 1 3 5 1 3 5 1 3 5 1 7 3
