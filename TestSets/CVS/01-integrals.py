from config import GEOMETRY_OPTIMIZED
from config import MOLE_ORB_TYPE
import os
from pyscf_util.misc.hforb_cacher import MoleHFOrbLoader
import numpy as np
from copy import deepcopy
import pyscf
from pyscf import tools

# argparse #

import argparse


def str2bool(v):
    """
    Convert string to boolean.

    Args:
        v: The input string.

    Returns:
        bool: The corresponding boolean value.
    """
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


parser = argparse.ArgumentParser(description="Dump iCIPT2 integrals")
parser.add_argument(
    "--basis", type=str, help="basis for the calculation", default="cc-pvtz"
)
parser.add_argument(
    "--x2c",
    type=str2bool,
    help="with X2C or not",
    default=True,
    nargs="?",
    const=True,
)
args = parser.parse_args()
print("Basis:", args.basis)
basis = args.basis
with_x2c = args.x2c

# load HF orb #

CVS_HF_loader = MoleHFOrbLoader(
    current_directory=os.path.dirname(os.path.realpath(__file__))
)

basis_dict = {
    "H": basis,
    "C": "unc-" + basis,
    "N": "unc-" + basis,
    "O": "unc-" + basis,
    "F": "unc-" + basis,
    "Si": "unc-" + basis,
    "P": "unc-" + basis,
    "S": "unc-" + basis,
    "Cl": "unc-" + basis,
}

FCIDUMP_FORMAT = "FCIDUMP_%s_%s_%s"  # mol_name, basis, task

MO_ENERGY_CUTOFF = 50.0

for mol_name in ["CO", "CO2", "CH4", "NH3", "H2O", "HF", "C2H2", "H2CO", "N2", "O2"]:

    if mol_name not in MOLE_ORB_TYPE.keys():
        continue

    mol, hf = CVS_HF_loader.get_mol_hf(
        mol_name, GEOMETRY_OPTIMIZED[mol_name], 0, 0, basis, with_x2c, True, False
    )
    # hf.analyze() #

    # extract mo_coeff to dump #

    mo_coeff = hf.mo_coeff
    mo_energy = hf.mo_energy

    idx = np.where(mo_energy < MO_ENERGY_CUTOFF)[0]
    mo_coeff = mo_coeff[:, idx].copy()

    print("mo_energy = ", mo_energy)
    print("mo_energy = ", mo_energy[idx])

    # reorder #

    norb_reorder = 0
    for key in MOLE_ORB_TYPE[mol_name]["orb_type"].keys():
        norb_reorder += len(MOLE_ORB_TYPE[mol_name]["orb_type"][key])

    print("norb_reorder:", norb_reorder)

    for task in MOLE_ORB_TYPE[mol_name]["task_type"]:
        print("dump task:", task)
        reorder = []
        reodering_indx = MOLE_ORB_TYPE[mol_name]["task_type"][task]["ordering"]
        print("reodering_indx:", reodering_indx)
        for key in reodering_indx:
            reorder += MOLE_ORB_TYPE[mol_name]["orb_type"][key]
        print("reorder:", reorder)

        mo_coeff_tmp = deepcopy(mo_coeff)
        mo_coeff_tmp[:, :norb_reorder] = mo_coeff[:, reorder]

        # dump #

        file_name = FCIDUMP_FORMAT % (mol_name, basis, task)
        if with_x2c:
            file_name = file_name + "_x2c"
        else:
            file_name = file_name + "_no_x2c"

        # transform integrals #

        int2e_full = pyscf.ao2mo.full(eri_or_mol=mol, mo_coeff=mo_coeff_tmp, aosym="s4")
        int1e_full = hf.get_hcore()
        int1e_full = np.dot(mo_coeff_tmp.T, np.dot(int1e_full, mo_coeff_tmp))
        energy_core = mol.enuc
        OrbSym = pyscf.symm.label_orb_symm(
            mol, mol.irrep_name, mol.symm_orb, mo_coeff_tmp
        )
        OrbSymID = [pyscf.symm.irrep_name2id(mol.groupname, x) for x in OrbSym]
        tools.fcidump.from_integrals(
            filename=file_name,
            h1e=int1e_full,
            h2e=int2e_full,
            nuc=energy_core,
            nmo=mo_coeff_tmp.shape[1],
            nelec=mol.nelectron,
            tol=1e-10,
            orbsym=OrbSymID,
        )
