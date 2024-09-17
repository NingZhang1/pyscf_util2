from config import MOLINFO
import os
from pyscf_util.misc.hforb_cacher import MoleHFOrbLoader
import pyscf
from pyscf import tools
import numpy as np
import argparse

mol_involved = list(MOLINFO.keys())
print("mol rerun for G2 test is ", mol_involved)

# argparse #

parser = argparse.ArgumentParser(description="Dump iCIPT2 integrals")
parser.add_argument(
    "--basis", type=str, help="basis for the calculation", default="cc-pvtz"
)
args = parser.parse_args()
print("Basis:", args.basis)
basis = args.basis

# load HF orb #

CVS_HF_loader = MoleHFOrbLoader(
    current_directory=os.path.dirname(os.path.realpath(__file__))
)

FCIDUMP_FORMAT = "FCIDUMP_%s_%s"  # mol_name, basis

for mol_name in mol_involved:

    mol, hf = CVS_HF_loader.get_mol_hf(
        mol_name,
        MOLINFO[mol_name]["geometry"],
        0,
        MOLINFO[mol_name]["spin"],
        basis,
        False,
        True,
        False,
    )
    hf.analyze()

    # dump integrals #

    file_name = FCIDUMP_FORMAT % (mol_name, basis)

    # transform integrals #

    mo_coeff = hf.mo_coeff

    int2e_full = pyscf.ao2mo.full(eri_or_mol=mol, mo_coeff=mo_coeff, aosym="s4")
    int1e_full = hf.get_hcore()
    int1e_full = np.dot(mo_coeff.T, np.dot(int1e_full, mo_coeff))
    energy_core = mol.enuc
    OrbSym = pyscf.symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, mo_coeff)
    OrbSymID = [pyscf.symm.irrep_name2id(mol.groupname, x) for x in OrbSym]
    tools.fcidump.from_integrals(
        filename=file_name,
        h1e=int1e_full,
        h2e=int2e_full,
        nuc=energy_core,
        nmo=mo_coeff.shape[1],
        nelec=mol.nelectron,
        tol=1e-10,
        orbsym=OrbSymID,
    )
