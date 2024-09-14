from config import MOLINFO
import os
from pyscf_util.misc.hforb_cacher import MoleHFOrbLoader
import pyscf
from pyscf import tools
import numpy as np
from pyscf_util.iCIPT2.iCIPT2 import kernel

mol_involved = list(MOLINFO.keys())
print("mol rerun for G2 test is ", mol_involved)

# argparse #

import argparse

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

CMIN = "5e-5 3e-5 2e-5 1e-5 7e-6 5e-6 3e-6 1.5e-6"

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

    # integrals #

    file_name = FCIDUMP_FORMAT % (mol_name, basis)

    nleft = MOLINFO[mol_name]["nleft"]
    nelec_val = MOLINFO[mol_name]["nelec"]

    kernel(
        True,
        task_name="%s_%s" % (basis, mol_name),
        fcidump=file_name,
        segment=MOLINFO[mol_name]["segment"] % (hf.mo_coeff.shape[1] - nleft),
        nelec_val=nelec_val,
        cmin=CMIN,
        Task=MOLINFO[mol_name]["task"],
        perturbation=1,
        rotatemo=1
    )
