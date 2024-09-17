from pyscf import gto, scf
from pyscf import tools
from pyscf_util.iCIPT2.iCIPT2_dooh import kernel
from pyscf_util.Integrals.integral_Dooh import FCIDUMP_Dooh

mol = gto.M(
    verbose=4,
    atom="""
            Cr   0.000000000000       0.000000000000      -0.84
            Cr   0.000000000000       0.000000000000       0.84
            """,
    basis="ccpvdz-dk",
    spin=0,
    charge=0,
    symmetry="dooh",
    unit="angstrom",
)
mol.build()
mf = scf.RHF(mol)
mf = scf.sfx2c(mf)
mf.kernel()
FCIDUMP_NAME = "FCIDUMP_Cr2_Dooh"

FCIDUMP_Dooh(mol, mf, FCIDUMP_NAME)

kernel(
    True,
    task_name="iCIPT2_Cr2_Dooh_CSF",
    fcidump=FCIDUMP_NAME,
    segment="10 8 6 6 %d 0" % (mol.nao - 30),
    nelec_val=12,
    cmin="1e-3 5e-4 1e-4 7e-5 5e-5 3e-5 1e-5",
    perturbation=1,
    Task="0 0 1 1",
)

kernel(
    False,
    task_name="iCIPT2_Cr2_Dooh_DET",
    fcidump=FCIDUMP_NAME,
    segment="10 8 6 6 %d 0" % (mol.nao - 30),
    nelec_val=12,
    cmin="1e-3 5e-4 1e-4 7e-5 5e-5 3e-5 1e-5",
    perturbation=1,
    Task="0 0 1 1",
)
