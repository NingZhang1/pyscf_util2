_CONFIG = {
    "Sc": {
        "sfx2c": True,
        0: {
            "mf": "mcscf",
            # description of symmetry #
            "spin": 1,
            "ici_state": [
                [1, 0, 2, [1, 1]],
                [1, 1, 1, [1]],
                [1, 2, 1, [1]],
                [1, 3, 1, [1]],
            ],
            "minimal_cas": {
                "norb": 6,
                "nelec": 3,
            },
            "fakescf": True,
            "fake_charge": 1,
            "fake_spin": 0,
            "fake_irrep_nelec": {
                "s+0": 8,
                "p+0": 4,
                "p-1": 4,
                "p+1": 4,
            },
            "cas_symm_d2h": {
                "ag": 3,
                "b1g": 1,
                "b2g": 1,
                "b3g": 1,
            },
        },
        1: {
            "mf": "rohf",
            "spin": 0,
        },
    },
    "Ti": {
        "sfx2c": True,
        0: {
            "mf": "mcscf",
            # description of symmetry #
            "spin": 2,
            "ici_state": [
                [2, 0, 1, [1]],
                [2, 1, 2, [1, 1]],
                [2, 2, 2, [1, 1]],
                [2, 3, 2, [1, 1]],
            ],
            "minimal_cas": {
                "norb": 6,
                "nelec": 4,
            },
            "fakescf": True,
            "fake_charge": 2,
            "fake_spin": 0,
            "fake_irrep_nelec": {
                "s+0": 8,
                "p+0": 4,
                "p-1": 4,
                "p+1": 4,
            },
            "cas_symm_d2h": {
                "ag": 3,
                "b1g": 1,
                "b2g": 1,
                "b3g": 1,
            },
        },
        1: {
            "mf": "mcscf",
            # description of symmetry #
            "spin": 3,
            "ici_state": [
                [3, 0, 1, [1]],
                [3, 1, 2, [1, 1]],
                [3, 2, 2, [1, 1]],
                [3, 3, 2, [1, 1]],
            ],
            "minimal_cas": {
                "norb": 6,
                "nelec": 3,
            },
            "fakescf": True,
            "fake_charge": 2,
            "fake_spin": 0,
            "fake_irrep_nelec": {
                "s+0": 8,
                "p+0": 4,
                "p-1": 4,
                "p+1": 4,
            },
            "cas_symm_d2h": {
                "ag": 3,
                "b1g": 1,
                "b2g": 1,
                "b3g": 1,
            },
        },
    },
    "V": {
        "sfx2c": True,
        0: {
            "mf": "mcscf",
            # description of symmetry #
            "spin": 3,
            "ici_state": [
                [3, 0, 1, [1]],
                [3, 1, 2, [1, 1]],
                [3, 2, 2, [1, 1]],
                [3, 3, 2, [1, 1]],
            ],
            "minimal_cas": {
                "norb": 6,
                "nelec": 5,
            },
            "fakescf": True,
            "fake_charge": 3,
            "fake_spin": 0,
            "fake_irrep_nelec": {
                "s+0": 8,
                "p+0": 4,
                "p-1": 4,
                "p+1": 4,
            },
            "cas_symm_d2h": {
                "ag": 3,
                "b1g": 1,
                "b2g": 1,
                "b3g": 1,
            },
        },
        1: {
            "mf": "mcscf",
            # description of symmetry #
            "spin": 4,
            "ici_state": [
                [4, 0, 1, [1]],
                [4, 1, 2, [1, 1]],
                [4, 2, 2, [1, 1]],
                [4, 3, 2, [1, 1]],
            ],
            "minimal_cas": {
                "norb": 6,
                "nelec": 4,
            },
            "fakescf": True,
            "fake_charge": 3,
            "fake_spin": 0,
            "fake_irrep_nelec": {
                "s+0": 8,
                "p+0": 4,
                "p-1": 4,
                "p+1": 4,
            },
            "cas_symm_d2h": {
                "ag": 3,
                "b1g": 1,
                "b2g": 1,
                "b3g": 1,
            },
        },
    },
    "Cr": {
        "sfx2c": True,
        0: {
            "mf": "rohf",
            "spin": 6,
        },
        1: {
            "mf": "rohf",
            "spin": 5,
        },
    },
    "Mn": {
        "sfx2c": True,
        0: {
            "mf": "rohf",
            "spin": 5,
        },
        1: {
            "mf": "rohf",
            "spin": 6,
        },
    },
    "Fe": {
        "sfx2c": True,
        0: {
            "mf": "mcscf",
            # description of symmetry #
            "spin": 4,
            "ici_state": [
                [4, 0, 2, [1, 1]],
                [4, 1, 1, [1]],
                [4, 2, 1, [1]],
                [4, 3, 1, [1]],
            ],
            "minimal_cas": {
                "norb": 6,
                "nelec": 8,
            },
            "fakescf": True,
            "fake_charge": 1,
            "fake_spin": 5,
            "fake_irrep_nelec": {
                "s+0": 8,
                "p+0": 4,
                "p-1": 4,
                "p+1": 4,
                "d+2": 1,
                "d+1": 1,
                "d+0": 1,
                "d-1": 1,
                "d-2": 1,
            },
            "cas_symm_d2h": {
                "ag": 3,
                "b1g": 1,
                "b2g": 1,
                "b3g": 1,
            },
        },
        1: {
            "mf": "mcscf",
            # description of symmetry #
            "spin": 5,
            "ici_state": [
                [5, 0, 2, [1, 1]],
                [5, 1, 1, [1]],
                [5, 2, 1, [1]],
                [5, 3, 1, [1]],
            ],
            "minimal_cas": {
                "norb": 6,
                "nelec": 7,
            },
            "fakescf": True,
            "fake_charge": 1,
            "fake_spin": 5,
            "fake_irrep_nelec": {
                "s+0": 8,
                "p+0": 4,
                "p-1": 4,
                "p+1": 4,
                "d+2": 1,
                "d+1": 1,
                "d+0": 1,
                "d-1": 1,
                "d-2": 1,
            },
            "cas_symm_d2h": {
                "ag": 3,
                "b1g": 1,
                "b2g": 1,
                "b3g": 1,
            },
        },
    },
    "Co": {
        "sfx2c": True,
        0: {
            "mf": "mcscf",
            # description of symmetry #
            "spin": 3,
            "ici_state": [
                [3, 0, 1, [1]],
                [3, 1, 2, [1, 1]],
                [3, 2, 2, [1, 1]],
                [3, 3, 2, [1, 1]],
            ],
            "minimal_cas": {
                "norb": 5,
                "nelec": 7,
            },
            "fakescf": True,
            "fake_charge": 2,
            "fake_spin": 5,
            "fake_irrep_nelec": {
                "s+0": 8,
                "p+0": 4,
                "p-1": 4,
                "p+1": 4,
                "d+2": 1,
                "d+1": 1,
                "d+0": 1,
                "d-1": 1,
                "d-2": 1,
            },
            "cas_symm_d2h": {
                "ag": 2,
                "b1g": 1,
                "b2g": 1,
                "b3g": 1,
            },
        },
        1: {
            "mf": "mcscf",
            # description of symmetry #
            "spin": 2,
            "ici_state": [
                [2, 0, 1, [1]],
                [2, 1, 2, [1, 1]],
                [2, 2, 2, [1, 1]],
                [2, 3, 2, [1, 1]],
            ],
            "minimal_cas": {
                "norb": 5,
                "nelec": 8,
            },
            "fakescf": True,
            "fake_charge": 4,
            "fake_spin": 5,
            "fake_irrep_nelec": {
                "s+0": 6,
                "p+0": 4,
                "p-1": 4,
                "p+1": 4,
                "d+2": 1,
                "d+1": 1,
                "d+0": 1,
                "d-1": 1,
                "d-2": 1,
            },
            "cas_symm_d2h": {
                "ag": 2,
                "b1g": 1,
                "b2g": 1,
                "b3g": 1,
            },
        },
    },
    "Ni": {
        "sfx2c": True,
        0: {
            "mf": "mcscf",
            # description of symmetry #
            "spin": 2,
            "ici_state": [
                [2, 0, 1, [1]],
                [2, 1, 2, [1, 1]],
                [2, 2, 2, [1, 1]],
                [2, 3, 2, [1, 1]],
            ],
            "minimal_cas": {
                "norb": 5,
                "nelec": 8,
            },
            "fakescf": True,
            "fake_charge": 3,
            "fake_spin": 5,
            "fake_irrep_nelec": {
                "s+0": 8,
                "p+0": 4,
                "p-1": 4,
                "p+1": 4,
                "d+2": 1,
                "d+1": 1,
                "d+0": 1,
                "d-1": 1,
                "d-2": 1,
            },
            "cas_symm_d2h": {
                "ag": 2,
                "b1g": 1,
                "b2g": 1,
                "b3g": 1,
            },
        },
        1: {
            "mf": "mcscf",
            # description of symmetry #
            "spin": 1,
            "ici_state": [
                [1, 0, 2, [1, 1]],
                [1, 1, 1, [1]],
                [1, 2, 1, [1]],
                [1, 3, 1, [1]],
            ],
            "minimal_cas": {
                "norb": 5,
                "nelec": 9,
            },
            "fakescf": True,
            "fake_charge": 0,
            "fake_spin": 0,
            "fake_irrep_nelec": {
                "s+0": 6,
                "p+0": 4,
                "p-1": 4,
                "p+1": 4,
                "d+2": 2,
                "d+1": 2,
                "d+0": 2,
                "d-1": 2,
                "d-2": 2,
            },
            "cas_symm_d2h": {
                "ag": 2,
                "b1g": 1,
                "b2g": 1,
                "b3g": 1,
            },
        },
    },
    "Cu": {
        "sfx2c": True,
        0: {
            "mf": "rohf",
            "spin": 1,
        },
        1: {
            "mf": "rohf",
            "spin": 0,
        },
    },
    "Zn": {
        "sfx2c": True,
        0: {
            "mf": "rohf",
            "spin": 0,
        },
        1: {
            "mf": "rohf",
            "spin": 1,
        },
    },
}
