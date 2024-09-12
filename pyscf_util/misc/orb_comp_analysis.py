from pyscf.data.elements import _std_symbol
from pyscf_util._atmMinCASOrb._orb_loader import (
    LoadAtmHFOrbAllInfo,
    LoadAtmHFOrbGivenType,
)
from pyscf_util._atmMinCASOrb._maingroup import _MO_CONFIG
import numpy
import numpy as np
from tabulate import tabulate


def Analysis_MO_Component(
    mol,
    mo_coeff,
    mo_energy,
    first_nmo=None,
    basis: dict[str, str] = None,
    with_sfx2c=None,
    verbose=True,
):
    # step 1 fetch each atm lable, which should be uniqie #
    atm_labels = []
    for atm_id in range(mol.natm):
        atm_symbol = mol.atom_symbol(atm_id)
        atm_labels.append(atm_symbol)
    assert len(atm_labels) == len(set(atm_labels))
    # step 2 determine ncomp #
    ncomp = 0
    comp_key = []
    for atm_label in atm_labels:
        ncomp += len(_MO_CONFIG[_std_symbol(atm_label)]["orb_type"])
        for orb_type in _MO_CONFIG[_std_symbol(atm_label)]["orb_type"]:
            key = "%4s %2s" % (atm_label, orb_type)
            comp_key.append(key)
    if verbose:
        print("ncomp = ", ncomp)
        print("comp_key = ", comp_key)
    # step 3 calculate the res #
    if first_nmo is None:
        first_nmo = mol.nao
    first_nmo = min(first_nmo, mol.nao)
    mo_coeff = mo_coeff[:, :first_nmo]
    res = numpy.ndarray((first_nmo, ncomp), dtype=numpy.float64)
    icomp = 0
    ovlp = mol.intor("int1e_ovlp")
    for atm_label in atm_labels:
        idx = mol.search_ao_label(["%s.*" % (atm_label)])
        assert idx[-1] + 1 - idx[0] == len(idx)
        for orb_type in _MO_CONFIG[_std_symbol(atm_label)]["orb_type"]:
            # fetch orb #
            orb = LoadAtmHFOrbGivenType(
                atm_label,
                0,
                basis[_std_symbol(atm_label)],
                with_sfx2c=with_sfx2c,
                orb_symbol=orb_type,
            )
            norb_loaded = orb.shape[1]
            orb_in_mol = numpy.zeros((mol.nao, norb_loaded), dtype=numpy.float64)
            orb_in_mol[idx[0] : idx[-1] + 1, :] = orb
            proj = numpy.dot(mo_coeff.T, numpy.dot(ovlp, orb_in_mol))
            proj = np.square(proj)
            row_sums = np.sum(proj, axis=1)
            res[:, icomp] = row_sums
            icomp += 1
    return {
        "nmo": first_nmo,
        "comp_key": comp_key,
        "comp": res,
        "mo_energy": mo_energy[:first_nmo],
    }


def print_dict_as_table(data):
    """
    Print the given dictionary as a formatted table, including MO energy.

    Parameters:
    data (dict): A dictionary containing 'nmo', 'comp_key', 'comp', and 'mo_energy' keys.
                 'comp' should be a 2D numpy array of shape (nmo, len(comp_key)).
                 'mo_energy' should be a 1D numpy array of length nmo.

    Returns:
    str: A string representation of the formatted table.
    """
    nmo = data["nmo"]
    comp_key = data["comp_key"]
    comp = data["comp"]
    mo_energy = data["mo_energy"]

    # Ensure comp is a 2D numpy array
    comp = np.array(comp).reshape(nmo, len(comp_key))

    # Ensure mo_energy is a 1D numpy array
    mo_energy = np.array(mo_energy).flatten()

    # Create headers
    headers = ["MO", "Energy"] + comp_key

    # Create table data
    table_data = []
    for i in range(nmo):
        row = [f"MO {i+1}", mo_energy[i]] + list(comp[i])
        table_data.append(row)

    # Generate the table
    table = tabulate(table_data, headers=headers, tablefmt="grid", floatfmt=".6f")

    return table
