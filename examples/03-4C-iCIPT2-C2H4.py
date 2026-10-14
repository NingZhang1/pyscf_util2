from pyscf import gto, scf

from pyscf_util.Relativisitc.integral_4C import (
    FCIDUMP_Rela4C,
)
from pyscf_util.Relativisitc.purify_spinor import purify_pes_v2, purify_pes_d2h

mol = gto.M(
    atom="""
C   0.000000   0.000000   0.669500
C   0.000000   0.000000  -0.669500
H   0.928900   0.000000   1.232100
H  -0.928900   0.000000   1.232100
H   0.928900   0.000000  -1.232100
H  -0.928900   0.000000  -1.232100
    """,
    basis="3-21g",
    verbose=5,
    charge=0,
    spin=0,
    symmetry="d2h",
)
mol.build()
mf = scf.dhf.RDHF(mol)
mf.conv_tol = 1e-12
mf.kernel()

FCIDUMP_NAME = "FCIDUMP_C2H4"

# 可以先尝试对其大分量 alpha 实部的对称性
# res = purify_pes_d2h(mol, mf, non_vanishing_tol=1e-12, max_iter=8)
# n2c = mol.nao_2c()
# mo_coeff = mf.mo_coeff.copy()
# mo_coeff[:, n2c:] = res
# mf.mo_coeff = mo_coeff

res = purify_pes_v2(mol, mf)
n2c = mol.nao_2c()
mo_coeff = mf.mo_coeff.copy()
mo_coeff[:, n2c:] = res
mf.mo_coeff = mo_coeff

FCIDUMP_Rela4C(
    mol,
    mf,
    False,
    filename=FCIDUMP_NAME,
    mode="outcore",
    debug=True,
    npes=28,
    IsComplex=False,
)
