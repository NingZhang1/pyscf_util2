from config import GEOMETRY_OPTIMIZED, CORE_CONFIG, SIGMA_CONFIG, PI_CONFIG
from config import NORB_ANALYSIS
import os
from pyscf_util.misc.hforb_cacher import MoleHFOrbLoader
from pyscf.data.elements import _std_symbol
from pyscf_util.misc.orb_comp_analysis import Analysis_MO_Component, print_dict_as_table

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
    current_directory=os.path.dirname(os.path.realpath(__file__)), _decontract_core=True
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

comp_dict = {}

for mol_name in ["CO", "CO2", "CH4", "NH3", "H2O", "HF", "C2H2", "H2CO", "N2", "O2"]:

    mol, hf = CVS_HF_loader.get_mol_hf(
        mol_name, GEOMETRY_OPTIMIZED[mol_name], 0, 0, basis, with_x2c, True, False
    )
    hf.analyze()

    # analysis #
    # determine nmo to analysis #
    first_nmo = 0
    for i in range(mol.natm):
        first_nmo += NORB_ANALYSIS[_std_symbol(mol.atom_symbol(i))]
    info = Analysis_MO_Component(
        mol, hf.mo_coeff, hf.mo_energy, first_nmo, basis_dict, with_sfx2c=True
    )
    comp_dict[mol_name] = info

for mol_name in [
    "SiH4",
    "PH3",
    "H2S",
    "HCl",
    "SO2",
    "H3COH",
    "Cl2",
    "NNO",
    "CH3CN",
    "HCN",
    "O3",
]:

    mol, hf = CVS_HF_loader.get_mol_hf(
        mol_name, GEOMETRY_OPTIMIZED[mol_name], 0, 0, basis, with_x2c, True, False
    )
    hf.analyze()

    # analysis #

    # determine nmo to analysis #
    first_nmo = 0
    for i in range(mol.natm):
        first_nmo += NORB_ANALYSIS[_std_symbol(mol.atom_symbol(i))]
    info = Analysis_MO_Component(
        mol, hf.mo_coeff, hf.mo_energy, first_nmo, basis_dict, with_sfx2c=True
    )
    # print(print_dict_as_table(info))
    comp_dict[mol_name] = info

# final print #

for key in comp_dict.keys():
    print("The orbital component of molecule %s" % (key))
    print(print_dict_as_table(comp_dict[key]))
    print("\n")
