_CONFIG = {
    "H": {
        0: {
            "mf": "rohf",
            "spin": 1,
        },
    },
    "C": {
        0: {
            "mf": "mcscf",
            # description of symmetry #
            "spin": 2,
            "state": [
                [2, "p+0", 1],
                [2, "p+1", 1],
                [2, "p-1", 1],
            ],
            "ici_state": [
                [2, 1, 1, [1]],
                [2, 2, 1, [1]],
                [2, 3, 1, [1]],
            ],
            "minimal_cas": {
                "norb": 3,
                "nelec": 2,
            },
            "cas_symm": {
                "p-1": 1,
                "p+0": 1,
                "p+1": 1,
            },
            # "core_symm": {
            #     "s+0": 2,
            # },
            "fakescf": True,
            "fake_charge": 2,
            "fake_spin": 0,
            "fake_irrep_nelec": {"s+0": 4},
            "cas_symm_d2h": {
                "b1u": 1,
                "b2u": 1,
                "b3u": 1,
            },
        },
    },
    "N": {
        0: {
            "mf": "rohf",
            "spin": 3,
        },
    },
    "O": {
        0: {
            "mf": "mcscf",
            # description of symmetry #
            "spin": 2,
            "state": [
                [2, "p+0", 1],
                [2, "p+1", 1],
                [2, "p-1", 1],
            ],
            "ici_state": [
                [2, 1, 1, [1]],
                [2, 2, 1, [1]],
                [2, 3, 1, [1]],
            ],
            "minimal_cas": {
                "norb": 3,
                "nelec": 4,
            },
            "cas_symm": {
                "p-1": 1,
                "p+0": 1,
                "p+1": 1,
            },
            "core_symm": {
                "s+0": 2,
            },
            "fakescf": True,
            "fake_charge": -2,
            "fake_spin": 0,
            "fake_irrep_nelec": {
                "s+0": 4,
                "p-1": 2,
                "p+0": 2,
                "p+1": 2,
            },
            "cas_symm_d2h": {
                "b1u": 1,
                "b2u": 1,
                "b3u": 1,
            },
        },
    },
    "F": {
        0: {
            "mf": "mcscf",
            # description of symmetry #
            "spin": 1,
            "state": [
                [1, "p+0", 1],
                [1, "p+1", 1],
                [1, "p-1", 1],
            ],
            "ici_state": [
                [1, 5, 1, [1]],
                [1, 6, 1, [1]],
                [1, 7, 1, [1]],
            ],
            "minimal_cas": {
                "norb": 3,
                "nelec": 5,
            },
            "cas_symm": {
                "p-1": 1,
                "p+0": 1,
                "p+1": 1,
            },
            "core_symm": {
                "s+0": 2,
            },
            "fakescf": True,
            "fake_charge": -1,
            "fake_spin": 0,
            "fake_irrep_nelec": {
                "s+0": 4,
                "p-1": 2,
                "p+0": 2,
                "p+1": 2,
            },
            "cas_symm_d2h": {
                "b1u": 1,
                "b2u": 1,
                "b3u": 1,
            },
        },
    },
    "Si": {
        0: {
            "mf": "mcscf",
            # description of symmetry #
            "spin": 2,
            "ici_state": [
                [2, 1, 1, [1]],
                [2, 2, 1, [1]],
                [2, 3, 1, [1]],
            ],
            "minimal_cas": {
                "norb": 3,
                "nelec": 2,
            },
            "cas_symm_d2h": {
                "b1u": 1,
                "b2u": 1,
                "b3u": 1,
            },
            "fakescf": True,
            "fake_charge": 2,
            "fake_spin": 0,
            "fake_irrep_nelec": {
                "s+0": 6,
                "p+0": 2,
                "p-1": 2,
                "p+1": 2,
            },
        },
    },
    "P": {
        0: {
            "mf": "rohf",
            "spin": 3,
        },
    },
    "S": {
        0: {
            "mf": "mcscf",
            # description of symmetry #
            "spin": 2,
            "ici_state": [
                [2, 1, 1, [1]],
                [2, 2, 1, [1]],
                [2, 3, 1, [1]],
            ],
            "minimal_cas": {
                "norb": 3,
                "nelec": 4,
            },
            "cas_symm_d2h": {
                "b1u": 1,
                "b2u": 1,
                "b3u": 1,
            },
            "fakescf": True,
            "fake_charge": -2,
            "fake_spin": 0,
            "fake_irrep_nelec": {
                "s+0": 6,
                "p+0": 4,
                "p-1": 4,
                "p+1": 4,
            },
        },
    },
    "Cl": {
        0: {
            "mf": "mcscf",
            # description of symmetry #
            "spin": 1,
            "state": [
                [1, "p+0", 1],
                [1, "p+1", 1],
                [1, "p-1", 1],
            ],
            "ici_state": [
                [1, 5, 1, [1]],
                [1, 6, 1, [1]],
                [1, 7, 1, [1]],
            ],
            "minimal_cas": {
                "norb": 3,
                "nelec": 5,
            },
            "fakescf": True,
            "fake_charge": -1,
            "fake_spin": 0,
            "fake_irrep_nelec": {
                "s+0": 6,
                "p-1": 4,
                "p+0": 4,
                "p+1": 4,
            },
            "cas_symm_d2h": {
                "b1u": 1,
                "b2u": 1,
                "b3u": 1,
            },
        },
    },
}

_MO_CONFIG = {
    "H": {
        "orb_type": ["1s", "2s", "2p"],
    },
    "C": {
        "orb_type": ["1s", "2s", "2p", "3s", "3p", "3d"],
    },
    "N": {
        "orb_type": ["1s", "2s", "2p", "3s", "3p", "3d"],
    },
    "O": {
        "orb_type": ["1s", "2s", "2p", "3s", "3p", "3d"],
    },
    "F": {
        "orb_type": ["1s", "2s", "2p", "3s", "3p", "3d"],
    },
    "Si": {
        "orb_type": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d"],
    },
    "P": {
        "orb_type": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d"],
    },
    "S": {
        "orb_type": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d"],
    },
    "Cl": {
        "orb_type": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d"],
    },
}
