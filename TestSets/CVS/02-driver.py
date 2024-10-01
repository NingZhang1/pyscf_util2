from config import GEOMETRY_OPTIMIZED
from config import MOLE_ORB_TYPE
import os
from pyscf_util.misc.hforb_cacher import MoleHFOrbLoader
import numpy as np
from copy import deepcopy
import pyscf
from pyscf import tools
from pyscf_util.iCIPT2.iCIPT2_CVS import kernel

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

    # extract mo_coeff to dump #

    mo_coeff = hf.mo_coeff
    mo_energy = hf.mo_energy

    idx = np.where(mo_energy < MO_ENERGY_CUTOFF)[0]
    mo_coeff = mo_coeff[:, idx].copy()

    print("mo_energy = ", mo_energy)
    print("mo_energy = ", mo_energy[idx])

    # reorder #

    for task in MOLE_ORB_TYPE[mol_name]["task_type"]:
        if task == "gt":
            continue
        print("run  task : ", task)
        file_name = FCIDUMP_FORMAT % (mol_name, basis, task)
        if with_x2c:
            file_name = file_name + "_x2c"
        else:
            file_name = file_name + "_no_x2c"

        print("iCIPT2 CVS")

        nleft = MOLE_ORB_TYPE[mol_name]["task_type"][task]["taskinfo"]["nleft"]
        nelec_val = MOLE_ORB_TYPE[mol_name]["task_type"][task]["taskinfo"]["nelec_val"]
        ici_task = MOLE_ORB_TYPE[mol_name]["task_type"][task]["taskinfo"]["task"]

        kernel(
            True,
            task_name="iCIPT2_CVS_%s_%s_%s" % (basis, mol_name, task),
            fcidump=file_name,
            segment=MOLE_ORB_TYPE[mol_name]["task_type"][task]["taskinfo"]["all"][
                "segment"
            ]
            % (mo_coeff.shape[1] - nleft),
            nelec_val=nelec_val,
            cmin="1e-4",
            Task=ici_task,
            perturbation=0,
        )

        # print("iCIPT2 CVS + CoreRelax")

        # kernel(
        #     True,
        #     task_name="iCIPT2_CVS_Relax_%s_%s_%s" % (basis, mol_name, task),
        #     fcidump=file_name,
        #     segment=MOLE_ORB_TYPE[mol_name]["task_type"][task]["taskinfo"]["all"][
        #         "segment"
        #     ]
        #     % (mo_coeff.shape[1] - nleft),
        #     nelec_val=nelec_val,
        #     cmin="1e-4",
        #     Task=ici_task,
        #     relaxcore=True,
        #     perturbation=0,
        # )
