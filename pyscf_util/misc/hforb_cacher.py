import os

# import numpy as np
from pyscf.data.elements import _std_symbol
from pyscf_util.misc.mole import get_mol
from pyscf_util.MeanField.scf import kernel as scf

# from pyscf_util.MeanField.mcscf import kernel as mcscf
from pyscf_util.MeanField.iciscf import kernel as iciscf
import re
import pickle
import copy

NAME_FORMAT = "%s_%d_%d_%s_%s.dat"  # molename,charge,spin,basis,x2c


class MoleHFOrbLoader:

    ORB_FORMAT_PATTERN = r"^\d+[spdfghi]$"

    @classmethod
    def _is_valid_orbformat(cls, orb_format):
        assert isinstance(orb_format, str)
        return bool(re.match(cls.ORB_FORMAT_PATTERN, orb_format))

    def __init__(self, current_directory, _decontract_core=False):
        if current_directory is None:
            self.current_directory = os.getcwd()
        else:
            self.current_directory = current_directory
        self.cache_dir_name = "_cached"
        # if not exist, create cache directory
        if not os.path.exists(
            os.path.join(self.current_directory, self.cache_dir_name)
        ):
            os.makedirs(os.path.join(self.current_directory, self.cache_dir_name))
        self._decontract_core = _decontract_core

    def get_mol(self, molename, geometry, charge, spin, basis, IsBohr=False):
        assert isinstance(molename, str)
        assert isinstance(charge, int)
        assert isinstance(basis, str)

        mol = get_mol(
            geometry,
            charge,
            spin,
            basis,
            symmetry=True,
            unit="Bohr" if IsBohr else "Angstrom",
        )
        if self._decontract_core:
            basis_new = {}
            for atm_id in range(mol.natm):
                atm_pure_symbol = _std_symbol(mol.atom_symbol(atm_id))
                if atm_pure_symbol != "H":
                    basis_new[mol.atom_symbol(atm_id)] = "unc-" + basis
                else:
                    basis_new[mol.atom_symbol(atm_id)] = basis
            mol = get_mol(
                geometry,
                charge,
                spin,
                basis_new,
                symmetry=True,
                unit="Bohr" if IsBohr else "Angstrom",
            )
            print("basis_new = ", basis_new)
            # exit(1)
        mol.build()
        return mol

    def load(
        self, molename, geometry, charge, spin, basis, x2c, IsBohr=False, rebuild=False
    ):
        assert isinstance(molename, str)
        assert isinstance(charge, int)
        assert isinstance(basis, str)

        if x2c:
            filename = NAME_FORMAT % (molename, charge, spin, basis, "X2C")
        else:
            filename = NAME_FORMAT % (molename, charge, spin, basis, "NOX2C")

        # if find file in cache directory
        if (
            os.path.exists(
                os.path.join(self.current_directory, self.cache_dir_name, filename)
            )
            and not rebuild
        ):
            return self._load_from_cache(filename)["mo_coeff"]
        else:
            info = self.build_info(molename, geometry, charge, spin, basis, x2c, IsBohr)
            return info["mo_coeff"]

    def load_all(
        self, molename, geometry, charge, spin, basis, x2c, IsBohr=False, rebuild=False
    ):
        assert isinstance(molename, str)
        assert isinstance(charge, int)
        assert isinstance(basis, str)

        if x2c:
            filename = NAME_FORMAT % (molename, charge, spin, basis, "X2C")
        else:
            filename = NAME_FORMAT % (molename, charge, spin, basis, "NOX2C")

        if (
            os.path.exists(
                os.path.join(self.current_directory, self.cache_dir_name, filename)
            )
            and not rebuild
        ):
            return self._load_from_cache(filename)
        else:
            return self.build_info(molename, geometry, charge, spin, basis, x2c, IsBohr)

    def get_mol_hf(
        self, molename, geometry, charge, spin, basis, x2c, IsBohr=False, rebuild=False
    ):
        mol = self.get_mol(molename, geometry, charge, spin, basis, IsBohr)
        hf = self._rhf(mol, with_sfx2c=x2c, run=False)
        info = self.load_all(
            molename, geometry, charge, spin, basis, x2c, IsBohr, rebuild
        )
        hf.mo_coeff = info["mo_coeff"]
        hf.mo_energy = info["mo_energy"]
        hf.mo_occ = info["mo_occ"]
        dm1 = hf.make_rdm1(mo_coeff=hf.mo_coeff, mo_occ=hf.mo_occ)
        hf.kernel(dm1)
        return mol, hf

    ### some helper func ###

    def _load_from_cache(self, filename):
        print(" ---- load from cache from %s ---- " % filename)
        with open(
            os.path.join(self.current_directory, self.cache_dir_name, filename), "rb"
        ) as f:
            return pickle.load(file=f)

    def _save_to_cache(self, filename, data):
        print(" ---- save to cache to %s ---- " % filename)
        data = copy.deepcopy(data)
        with open(
            os.path.join(self.current_directory, self.cache_dir_name, filename), "wb"
        ) as f:
            pickle.dump(obj=data, file=f)
        return None

    def build_info(
        self, molename, geometry, charge, spin, basis, x2c=False, IsBohr=False
    ):

        with_sfx2c = x2c

        # hf info #

        # mol = self._get_mol(geometry, charge, spin, basis, IsBohr)
        mol = self.get_mol(molename, geometry, charge, spin, basis, IsBohr)
        hf = self._rhf(mol, with_sfx2c=with_sfx2c)
        info = {
            "ovlp": hf.get_ovlp(),
            # "mol": mol,
            "mo_coeff": hf.mo_coeff,
            "mo_energy": hf.mo_energy,
            "mo_occ": hf.mo_occ,
        }
        print("mo_energy", info["mo_energy"])
        # filename = NAME_FORMAT % (atm, charge, basis)
        if with_sfx2c:
            filename = NAME_FORMAT % (molename, charge, spin, basis, "X2C")
        else:
            filename = NAME_FORMAT % (molename, charge, spin, basis, "NOX2C")
        self._save_to_cache(filename, info)
        return info

    # def _get_mol(self, geometry, charge, spin, basis, IsBohr=False):
    #     xyz = geometry
    #     if not IsBohr:
    #         mol = get_mol(xyz, charge, spin, basis, symmetry=True)
    #     else:
    #         mol = get_mol(xyz, charge, spin, basis, symmetry=True, unit="Bohr")
    #     mol.build()
    #     return mol

    def _rhf(self, mol, with_sfx2c=False, irrep_nelec=None, run=True):
        return scf(mol, sfx1e=with_sfx2c, irrep_nelec=irrep_nelec, run=run)
