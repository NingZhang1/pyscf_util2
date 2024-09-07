import pyscf


def kernel(mol, irrep_nelec=None, sfx1e=False, newton=False, run=True):
    """Run SCF calculation given a molecule object
    Args:
        mol: a molecule object
        sfx1e: whether to use sfx2c
        newton: whether to use newton solver

    Kwargs:

    Returns:
        my_hf: a (runned) pyscf SCF object
    """

    my_hf = pyscf.scf.ROHF(mol)
    if sfx1e:
        my_hf = pyscf.scf.sfx2c(my_hf)
    if newton:
        my_hf = pyscf.scf.newton(my_hf)
    if irrep_nelec is not None:
        my_hf.irrep_nelec = irrep_nelec
    if run:
        my_hf.kernel()
    else:
        my_hf.build()
    return my_hf


scf = kernel
