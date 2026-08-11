from pyscf import gto, scf
from pyscf import tools
from pyscf_util.iCIPT2.iCIPT2 import kernel

mol = gto.M(
    verbose=4,
    atom="""
            H   0.000000000000       1.000000000000      0.00000
            H   1.000000000000       0.000000000000      0.00000
            O   0.000000000000       0.000000000000      0.00000
            """,
    basis="cc-pvdz",
    spin=0,
    charge=0,
    symmetry="c2v",
)
mol.build()
mf = scf.RHF(mol)
mf.kernel()

mf.analyze()

mf.MP2().kernel()

FCIDUMP_NAME = "H2O.FCIDUMP"

tools.fcidump.from_scf(mf, FCIDUMP_NAME, 1e-10)


