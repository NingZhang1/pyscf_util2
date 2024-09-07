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

NAME_FORMAT = "%s_%d_%s.dat"  # atm,charge,basis

### atm MINCAS config ###

# TODO: change to a json file

from pyscf_util._atmMinCASOrb._maingroup import _CONFIG as _CONFIG


class AtmHFOrbLoader:

    ORB_FORMAT_PATTERN = r"^\d+[spdfghi]$"

    @classmethod
    def _is_valid_orbformat(cls, orb_format):
        assert isinstance(orb_format, str)
        return bool(re.match(cls.ORB_FORMAT_PATTERN, orb_format))

    def __init__(self):
        self.current_directory = os.path.dirname(os.path.realpath(__file__))
        self.cache_dir_name = "_cached"
        # if not exist, create cache directory
        if not os.path.exists(
            os.path.join(self.current_directory, self.cache_dir_name)
        ):
            os.makedirs(os.path.join(self.current_directory, self.cache_dir_name))

    def load(self, atm, charge, basis):
        assert isinstance(atm, str)
        atm = _std_symbol(atm)
        assert isinstance(charge, int)
        assert isinstance(basis, str)

        filename = NAME_FORMAT % (atm, charge, basis)

        # if find file in cache directory
        if os.path.exists(
            os.path.join(self.current_directory, self.cache_dir_name, filename)
        ):
            return self._load_from_cache(filename)["mo_coeff"]
        else:
            info = self.build_info(atm, charge, basis)
            return info["mo_coeff"]

    def load_all(self, atm, charge, basis):
        assert isinstance(atm, str)
        atm = _std_symbol(atm)
        assert isinstance(charge, int)
        assert isinstance(basis, str)
        filename = NAME_FORMAT % (atm, charge, basis)
        if os.path.exists(
            os.path.join(self.current_directory, self.cache_dir_name, filename)
        ):
            return self._load_from_cache(filename)
        else:
            return self.build_info(atm, charge, basis)

    # def load_orb(self, atm, charge, basis, orb_symbol: list[str], cutoff=0.2):
    #     assert isinstance(atm, str)
    #     atm = _std_symbol(atm)
    #     assert isinstance(charge, int)
    #     assert isinstance(basis, str)
    #     assert isinstance(orb_symbol, list)
    #     assert all([self._is_valid_orbformat(x) for x in orb_symbol])
    #     res = {x: None for x in orb_symbol}
    #     filename = NAME_FORMAT % (atm, charge, basis)
    #     return

    ### some helper func ###

    def _load_from_cache(self, filename):
        print(" ---- load from cache from %s ---- " % filename)
        with open(
            os.path.join(self.current_directory, self.cache_dir_name, filename), "rb"
        ) as f:
            return pickle.load(f)

    def _save_to_cache(self, filename, data):
        print(" ---- save to cache to %s ---- " % filename)
        with open(
            os.path.join(self.current_directory, self.cache_dir_name, filename), "wb"
        ) as f:
            pickle.dump(data, f)
        return None

    def build_info(self, atm, charge, basis, with_sfx2c=None):

        if with_sfx2c is None:
            if "sfx2c" in _CONFIG[atm]:
                with_sfx2c = _CONFIG[atm]["sfx2c"]
            else:
                with_sfx2c = False

        # hf info #

        if _CONFIG[atm][charge]["mf"] == "rohf":
            mol = self._get_mol(atm, charge, _CONFIG[atm][charge]["spin"], basis)
            hf = self._rhf(mol, with_sfx2c=with_sfx2c)
            info = {
                "ovlp": hf.get_ovlp(),
                "mol": mol,
                "mo_coeff": hf.mo_coeff,
                "mo_energy": hf.mo_energy,
            }
        else:
            assert (
                _CONFIG[atm][charge]["mf"] == "mcscf"
            )  # TODO: it seems to be problematic when doing pyscf with high sym
            mol = self._get_mol(atm, charge, _CONFIG[atm][charge]["spin"], basis)
            if _CONFIG[atm][charge]["fakescf"]:
                mol_fake = self._get_mol(
                    atm,
                    _CONFIG[atm][charge]["fake_charge"],
                    _CONFIG[atm][charge]["fake_spin"],
                    basis,
                )
                hf_fake = self._rhf(
                    mol_fake,
                    irrep_nelec=_CONFIG[atm][charge]["fake_irrep_nelec"],
                    with_sfx2c=with_sfx2c,
                )
                init_guess = copy.deepcopy(hf_fake.mo_coeff)
                mf = self._rhf(mol, run=False)
            else:
                mf = self._rhf(mol)
                init_guess = copy.deepcopy(mf.mo_coeff)
            # perform mcscf #
            mol.symmetry = "D2h"
            mol.build()
            mf = self._rhf(mol, run=False)
            mc = self._mcscf(
                mol,
                mf,
                _CONFIG[atm][charge]["minimal_cas"]["norb"],
                _CONFIG[atm][charge]["minimal_cas"]["nelec"],
                _CONFIG[atm][charge]["ici_state"],
                init_guess=init_guess,
                cas_list=_CONFIG[atm][charge]["cas_symm_d2h"],
            )
            info = {
                "ovlp": mf.get_ovlp(),
                "mol": mol,
                "mo_coeff": mc.mo_coeff,
                "mo_energy": mc.mo_energy,
            }
        print("mo_energy", info["mo_energy"])
        filename = NAME_FORMAT % (atm, charge, basis)
        self._save_to_cache(filename, info)
        return info

    def _get_mol(self, atm, charge, spin, basis):
        xyz = "%s 0 0 0" % atm
        mol = get_mol(xyz, charge, spin, basis, symmetry=True)
        mol.build()
        return mol

    def _rhf(self, mol, with_sfx2c=False, irrep_nelec=None, run=True):

        return scf(mol, sfx1e=with_sfx2c, irrep_nelec=irrep_nelec, run=run)

    def _mcscf(self, mol, mf, norb, nelec, state, init_guess=None, cas_list=None):

        return iciscf(
            mol,
            mf,
            nelec,
            norb,
            _ici_state=state,
            _mo_init=init_guess,
            _cas_list=cas_list,
        )


LOADER = AtmHFOrbLoader()


def LoadAtmHFOrb(atm, charge, basis):
    return LOADER.load(atm, charge, basis)


def LoadAtmHFOrbAllInfo(atm, charge, basis):
    return LOADER.load_all(atm, charge, basis)


if __name__ == "__main__":
    print(LoadAtmHFOrb("C", 0, "cc-pvdz"))
    print(LoadAtmHFOrb("C", 0, "cc-pvdz"))
    print(LoadAtmHFOrb("N", 0, "cc-pvdz"))
    print(LoadAtmHFOrb("N", 0, "cc-pvdz"))
    print(LoadAtmHFOrb("O", 0, "cc-pvdz"))
    print(LoadAtmHFOrb("O", 0, "cc-pvdz"))
    print(LoadAtmHFOrb("F", 0, "cc-pvdz"))
    print(LoadAtmHFOrb("F", 0, "cc-pvdz"))
    print(LoadAtmHFOrb("Si", 0, "cc-pvdz"))
    print(LoadAtmHFOrb("P", 0, "cc-pvdz"))
    print(LoadAtmHFOrb("S", 0, "cc-pvdz"))
    print(LoadAtmHFOrb("Cl", 0, "cc-pvdz"))
    print(LoadAtmHFOrbAllInfo("F", 0, "cc-pvdz"))
