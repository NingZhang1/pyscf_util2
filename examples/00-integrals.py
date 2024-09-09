from pyscf_util._atmMinCASOrb._orb_loader import LoadAtmHFOrb
from pyscf_util.MeanField.scf import kernel as scf
from pyscf_util.MeanField.mcscf import kernel as mcscf
from pyscf_util.misc.mole import get_mol
from pyscf_util._atmMinCASOrb._transitionmetal import _CONFIG as _TM_CONFIG
from pyscf_util.Integrals.integral_CASCI import dump_heff_casci

# use argparse to parse command line arguments basis #

import argparse

parser = argparse.ArgumentParser(description="Dump CASCI integrals")
parser.add_argument(
    "--basis", type=str, help="basis for the calculation", default="cc-pvtz"
)
args = parser.parse_args()
print("Basis:", args.basis)

# start the calculation #

FCIDUMP_FORMAT = "FCIDUMP_%s_%d_%s"

for atm in ["Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"]:
    for charge in [0, 1]:
        mo_coeff = LoadAtmHFOrb(atm, charge, args.basis)
        xyz = "%s 0 0 0" % atm
        mol = get_mol(
            xyz, charge, _TM_CONFIG[atm][charge]["spin"], args.basis, symmetry="D2h"
        )
        # build hf obj #
        scf_obj = scf(mol, sfx1e=True, run=False)
        scf_obj.mo_coeff = mo_coeff
        # build mcscf obj #
        try:
            mcscf_obj, _ = mcscf(
                mol,
                scf_obj,
                _TM_CONFIG[atm][charge]["minimal_cas"]["nelec"],
                _TM_CONFIG[atm][charge]["minimal_cas"]["norb"],
                _mo_init=mo_coeff,
                _run_mcscf=False,
            )
        except:
            mcscf_obj, _ = mcscf(
                mol,
                scf_obj,
                0,
                0,
                _mo_init=mo_coeff,
                _run_mcscf=False,
            )
        # dump CASCI ints #
        dump_heff_casci(
            mol,
            mcscf_obj,
            mo_coeff[:, :5],
            mo_coeff[:, 5:],
            FCIDUMP_FORMAT % (atm, charge, basis),
        )
