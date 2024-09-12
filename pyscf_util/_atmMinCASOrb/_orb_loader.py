import os
import numpy as np

from pyscf.data.elements import _std_symbol
from pyscf_util.misc.mole import get_mol
from pyscf_util.MeanField.scf import kernel as scf

# from pyscf_util.MeanField.mcscf import kernel as mcscf
from pyscf_util.MeanField.iciscf import kernel as iciscf
import re
import pickle
import copy

NAME_FORMAT = "%s_%d_%s"  # atm,charge,basis

### atm MINCAS config ###

# TODO: change to a json file

from pyscf_util._atmMinCASOrb._maingroup import _CONFIG as _CONFIG
from pyscf_util._atmMinCASOrb._transitionmetal import _CONFIG as _TM_CONFIG
from pyscf_util._atmMinCASOrb._maingroup import _MO_CONFIG

_CONFIG.update(_TM_CONFIG)

##########################
######## HELPER ##########


def classify_orbitals(mo_energy, threshold=1e-5, verbose=False):

    # 初始化分组
    groups = []
    current_group = [mo_energy[0]]

    # 根据阈值对能量值进行分组
    for energy in mo_energy[1:]:
        if energy - current_group[-1] <= threshold:
            current_group.append(energy)
        else:
            groups.append(current_group)
            current_group = [energy]
    groups.append(current_group)

    principle_number = {"s": 1, "p": 2, "d": 3, "f": 4, "g": 5, "h": 6, "i": 7}

    # 确定每组的轨道类型
    orbital_types = []
    for group in groups:
        count = len(group)
        if count == 1:
            orbital_types.extend([f'{principle_number["s"]}s'] * count)
            principle_number["s"] += 1
        elif count == 3:
            orbital_types.extend([f'{principle_number["p"]}p'] * count)
            principle_number["p"] += 1
        elif count == 5:
            orbital_types.extend([f'{principle_number["d"]}d'] * count)
            principle_number["d"] += 1
        elif count == 7:
            orbital_types.extend([f'{principle_number["f"]}f'] * count)
            principle_number["f"] += 1
        elif count == 9:
            orbital_types.extend([f'{principle_number["g"]}g'] * count)
            principle_number["g"] += 1
        elif count == 11:
            orbital_types.extend([f'{principle_number["h"]}h'] * count)
            principle_number["h"] += 1
        elif count == 13:
            orbital_types.extend([f'{principle_number["i"]}i'] * count)
            principle_number["i"] += 1
        else:
            orbital_types.extend(["unknown"] * count)

    # 打印结果
    if verbose:
        print("能量\t\t轨道类型")
        for energy, orbital_type in zip(mo_energy, orbital_types):
            print(f"{energy:15.6f}\t\t{orbital_type}")

    return orbital_types


##########################


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
        self._cached = {}

    def load(self, atm, charge, basis, with_sfx2c=None):
        assert isinstance(atm, str)
        atm = _std_symbol(atm)
        assert isinstance(charge, int)
        assert isinstance(basis, str)

        filename = self._get_filename(atm, charge, basis, with_sfx2c)

        # if find file in cache directory
        if os.path.exists(
            os.path.join(self.current_directory, self.cache_dir_name, filename)
        ):
            return self._load_from_cache(filename)["mo_coeff"]
        else:
            info = self.build_info(atm, charge, basis, with_sfx2c)
            return info["mo_coeff"]

    def load_all(self, atm, charge, basis, with_sfx2c=None):
        assert isinstance(atm, str)
        atm = _std_symbol(atm)
        assert isinstance(charge, int)
        assert isinstance(basis, str)

        filename = self._get_filename(atm, charge, basis, with_sfx2c)

        if os.path.exists(
            os.path.join(self.current_directory, self.cache_dir_name, filename)
        ):
            return self._load_from_cache(filename)
        else:
            return self.build_info(atm, charge, basis, with_sfx2c)

    def load_orb(
        self, atm, charge, basis, with_sfx2c=None, orb_symbol: list[str] = None
    ):
        assert isinstance(atm, str)
        atm = _std_symbol(atm)
        assert isinstance(charge, int)
        assert isinstance(basis, str)
        assert isinstance(orb_symbol, str) or isinstance(orb_symbol, list)
        assert orb_symbol is not None
        if isinstance(orb_symbol, str):
            orb_symbol = [orb_symbol]
        assert all(self._is_valid_orbformat(orb) for orb in orb_symbol)
        assert all([x in _MO_CONFIG[_std_symbol(atm)]["orb_type"] for x in orb_symbol])

        # print("with_sfx2c = ", with_sfx2c)
        info = self.load_all(atm, charge, basis, with_sfx2c)
        orbs_type = classify_orbitals(info["mo_energy"])
        indx = []
        for orb_type in orb_symbol:
            tmp = np.where(np.array(orbs_type) == orb_type)[0]
            indx.extend(tmp)
        return info["mo_coeff"][:, indx]

    ### some helper func ###

    def _get_filename(self, atm, charge, basis, with_sfx2c=None):
        if with_sfx2c is None:
            if "sfx2c" in _CONFIG[atm]:
                with_sfx2c = _CONFIG[atm]["sfx2c"]
            else:
                with_sfx2c = False

        filename = NAME_FORMAT % (atm, charge, basis)
        if with_sfx2c:
            filename += "_X2C.dat"
        else:
            filename += "_NOX2C.dat"
        return filename

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

    def build_info(self, atm, charge, basis, with_sfx2c=None):

        if with_sfx2c is None:
            if "sfx2c" in _CONFIG[atm]:
                with_sfx2c = _CONFIG[atm]["sfx2c"]
            else:
                with_sfx2c = False

        # hf info #

        if _CONFIG[atm][charge]["mf"] == "rohf":
            mol = self._get_mol(atm, charge, _CONFIG[atm][charge]["spin"], basis)
            if "irrep_nelec" in _CONFIG[atm][charge].keys():
                # print("irrep_nelec", _CONFIG[atm][charge]["irrep_nelec"])
                hf = self._rhf(
                    mol,
                    irrep_nelec=_CONFIG[atm][charge]["irrep_nelec"],
                    with_sfx2c=with_sfx2c,
                )
            else:
                hf = self._rhf(mol, with_sfx2c=with_sfx2c)
            info = {
                "ovlp": hf.get_ovlp(),
                # "mol": mol,
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
                # print(hf_fake.mo_energy)
                # exit(1)
                init_guess = copy.deepcopy(hf_fake.mo_coeff)
                mf = self._rhf(mol, run=False, with_sfx2c=with_sfx2c)
            else:
                mf = self._rhf(mol, with_sfx2c=with_sfx2c)
                init_guess = copy.deepcopy(mf.mo_coeff)
            # perform mcscf #
            mol.symmetry = "D2h"
            mol.build()
            mf = self._rhf(mol, run=False, with_sfx2c=with_sfx2c)
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
                # "mol": mol,
                "mo_coeff": mc.mo_coeff,
                "mo_energy": mc.mo_energy,
            }
        print("mo_energy", info["mo_energy"])
        filename = self._get_filename(atm, charge, basis, with_sfx2c)
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
            _mc_max_macro=20,
        )


LOADER = AtmHFOrbLoader()


def LoadAtmHFOrb(atm, charge, basis, with_sfx2c=None, rerun=False):
    if not rerun:
        return LOADER.load(atm, charge, basis, with_sfx2c)
    else:
        return LOADER.build_info(atm, charge, basis, with_sfx2c)


def LoadAtmHFOrbAllInfo(atm, charge, basis, with_sfx2c=None, rerun=False):
    if not rerun:
        return LOADER.load_all(atm, charge, basis, with_sfx2c)
    else:
        return LOADER.build_info(atm, charge, basis, with_sfx2c)


def LoadAtmHFOrbGivenType(
    atm, charge, basis, with_sfx2c=None, orb_symbol=None, rerun=False
):
    if not rerun:
        return LOADER.load_orb(atm, charge, basis, with_sfx2c, orb_symbol)
    else:
        LOADER.build_info(atm, charge, basis, with_sfx2c)
        return LOADER.load_orb(atm, charge, basis, with_sfx2c, orb_symbol)


if __name__ == "__main__":
    # print(LoadAtmHFOrb("C", 0, "cc-pvdz"))
    # print(LoadAtmHFOrb("C", 0, "cc-pvdz"))
    # print(LoadAtmHFOrb("N", 0, "cc-pvdz"))
    # print(LoadAtmHFOrb("N", 0, "cc-pvdz"))
    # print(LoadAtmHFOrb("O", 0, "cc-pvdz"))
    # print(LoadAtmHFOrb("O", 0, "cc-pvdz"))
    # print(LoadAtmHFOrb("F", 0, "cc-pvdz"))
    # print(LoadAtmHFOrb("F", 0, "cc-pvdz"))
    # print(LoadAtmHFOrb("Si", 0, "cc-pvdz"))
    # print(LoadAtmHFOrb("P", 0, "cc-pvdz"))
    # print(LoadAtmHFOrb("S", 0, "cc-pvdz"))
    # print(LoadAtmHFOrb("Cl", 0, "cc-pvdz"))
    # print(LoadAtmHFOrbAllInfo("F", 0, "cc-pvdz"))
    # print(LoadAtmHFOrb("Sc", 0, "cc-pvtz"))
    # print(LoadAtmHFOrb("Sc", 1, "cc-pvtz"))
    # print(LoadAtmHFOrb("Sc", 1, "cc-pvtz", rerun=True))
    # # print(LoadAtmHFOrb("Ti", 0, "cc-pvtz", rerun=True))
    # print(LoadAtmHFOrb("Ti", 1, "cc-pvtz", rerun=True))
    # print(LoadAtmHFOrb("V", 1, "cc-pvtz", rerun=True))
    # print(LoadAtmHFOrb("Co", 0, "cc-pvtz", rerun=True))
    # print(LoadAtmHFOrb("Co", 1, "cc-pvtz", rerun=True))
    # print(LoadAtmHFOrb("Ni", 1, "cc-pvdz", rerun=True))
    # print(LoadAtmHFOrb("Fe", 0, "cc-pvtz", rerun=True))
    # print(LoadAtmHFOrb("Fe", 1, "cc-pvtz", rerun=True))
    # print(LoadAtmHFOrb("Co", 0, "cc-pvtz", rerun=True))
    # print(LoadAtmHFOrb("Co", 1, "cc-pvtz", rerun=True))
    # print(LoadAtmHFOrb("Ni", 0, "cc-pvtz", rerun=True))
    # print(LoadAtmHFOrb("Ni", 1, "cc-pvtz", rerun=True))
    print(LoadAtmHFOrb("Cr", 1, "cc-pvtz", rerun=True))
