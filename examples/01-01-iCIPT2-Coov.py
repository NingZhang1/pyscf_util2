from pyscf import gto, scf
from pyscf import tools
from pyscf_util.iCIPT2.iCIPT2_coov import kernel
from pyscf_util.Integrals.integral_Coov import (
    FCIDUMP_Coov,
)

mol = gto.M(
    verbose=0,
    atom="""
            C   0.000000000000       0.000000000000      -0.621265
            O   0.000000000000       0.000000000000       0.621265
            """,
    basis={"C": "cc-pvdz", "O": "cc-pvdz"},
    spin=0,
    charge=0,
    symmetry="coov",
)
mol.build()
mf = scf.RHF(mol)
mf.kernel()

FCIDUMP_NAME = "FCIDUMP_CO_COOV"

FCIDUMP_Coov(mol, mf, FCIDUMP_NAME)

kernel(
    True,
    task_name="iCIPT2_CO_Coov",
    fcidump=FCIDUMP_NAME,
    segment="2 2 4 4 %d 0" % (mol.nao - 12),
    nelec_val=6,
    cmin="1e-3",
    perturbation=1,
    Task="0 0 1 1",
    end_with=".PrimeSpace",
)

kernel(
    True,
    task_name="iCIPT2_CO_Coov",
    fcidump=FCIDUMP_NAME,
    segment="2 2 4 4 %d 0" % (mol.nao - 12),
    nelec_val=6,
    cmin="1e-4",
    perturbation=1,
    Task="0 0 1 1",
)
