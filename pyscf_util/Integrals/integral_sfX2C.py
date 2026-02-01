import pyscf
import numpy
from functools import reduce
from pyscf import tools
import tempfile
import h5py
from pyscf.ao2mo import outcore
from pyscf.mcscf.casci import get_fock
from pyscf_util.misc.misc import _combine4, _combine2
from pyscf_util.File import file_cmoao

from functools import reduce
import numpy
from pyscf.tools import fcidump


def fcidump_sfx2c(mol, scf, mo_coeff, filename="FCIDUMP", tol=1e-10):

    nelec = mol.nelectron
    ms = 0
    # tol = 1e-10
    nuc = mol.get_enuc()

    nao = mo_coeff.shape[1]

    h1e = reduce(numpy.dot, (mo_coeff.T, scf.get_hcore(), mo_coeff))
    h1e = h1e[:nao, :nao]

    # print(h1e)

    int2e_full = pyscf.ao2mo.full(eri_or_mol=mol, mo_coeff=mo_coeff[:, :nao], aosym="4")
    # int2e_full = pyscf.ao2mo.restore(8, int2e_full.copy(), nao)

    OrbSym = pyscf.symm.label_orb_symm(
        mol, mol.irrep_name, mol.symm_orb, mo_coeff[:, :nao]
    )
    OrbSymID = [pyscf.symm.irrep_name2id(mol.groupname, x) for x in OrbSym]

    fcidump.from_integrals(
        filename, h1e, int2e_full, nao, nelec, nuc, ms, OrbSymID, tol
    )
