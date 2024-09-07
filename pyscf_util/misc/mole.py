import pyscf


def get_orbsym(mol, mocoeff):

    OrbSym = pyscf.symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, mocoeff)
    OrbSymID = [pyscf.symm.irrep_name2id(mol.groupname, x) for x in OrbSym]

    return OrbSymID, OrbSym


def get_mol(
    xyz, charge=0, spin=0, basis="6-31G(d)", symmetry="", verbose=4, unit="angstorm"
):
    mol = pyscf.gto.M(
        verbose=verbose,
        atom=xyz,
        basis=basis,
        spin=spin,
        charge=charge,
        symmetry=symmetry,
        unit=unit,
    )
    mol.build()
    return mol
