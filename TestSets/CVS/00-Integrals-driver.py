from config import GEOMETRY_OPTIMIZED, CORE_CONFIG, SIGMA_CONFIG, PI_CONFIG
import os
from pyscf_util.misc.hforb_cacher import MoleHFOrbLoader
from pyscf.data.elements import _std_symbol

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

for mol_name in ["CO", "CO2", "CH4", "NH3", "H2O", "HF", "C2H2", "H2CO", "N2", "O2"]:

    mol, hf = CVS_HF_loader.get_mol_hf(
        mol_name, GEOMETRY_OPTIMIZED[mol_name], 0, 0, basis, with_x2c, True, False
    )
    hf.analyze()

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
