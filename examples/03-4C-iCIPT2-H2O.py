from pyscf import gto, scf

from pyscf_util.Relativisitc.integral_4C import (
    FCIDUMP_Rela4C,
)
from pyscf_util.Relativisitc.purify_spinor import purify_pes_v2, purify_pes_d2h

mol = gto.M(
    atom="""
        O   0.000000   0.000000   0.000000
        H   0.000000  -0.757000   0.587000
        H   0.000000   0.757000   0.587000
    """,
    basis="3-21g",
    verbose=5,
    charge=0,
    spin=0,
    symmetry="c2v",
)
mol.build()
mf = scf.dhf.RDHF(mol)
mf.conv_tol = 1e-12
mf.kernel()

FCIDUMP_NAME = "FCIDUMP_H2O"

res = purify_pes_v2(mol, mf)  # 纯化，但不保证正确
n2c = mol.nao_2c()
mo_coeff = mf.mo_coeff.copy()
mo_coeff[:, n2c:] = res
mf.mo_coeff = mo_coeff

coulomb, _ = FCIDUMP_Rela4C(
    mol,
    mf,
    False,
    filename=FCIDUMP_NAME,
    mode="incore",
    debug=True,
    npes=20,
    IsComplex=False,
)

FCIDUMP_Rela4C(
    mol,
    mf,
    False,
    filename=FCIDUMP_NAME + "_nosymm",
    mode="incore",  # 太大的话可以 "outcore"
    debug=True,  # 也可以为 True
    npes=20,  # 限制有多少个 pes spinor (kramers pairs x 2)
    IsComplex=True,
)  # 产生不用对称性做 4c-icipt2 的积分

for i in range(20):
    for j in range(20):
        for k in range(20):
            for l in range(20):
                n_odd = 0
                if i % 2 == 1:
                    n_odd += 1
                if j % 2 == 1:
                    n_odd += 1
                if k % 2 == 1:
                    n_odd += 1
                if l % 2 == 1:
                    n_odd += 1
                if n_odd % 2 == 1:
                    if abs(coulomb[i, j, k, l]) >= 1e-7:
                        print(
                            "warning : coulomb[%3d,%3d,%3d,%3d] = %20.16e not zero"
                            % (i, j, k, l, coulomb[i, j, k, l])
                        )
