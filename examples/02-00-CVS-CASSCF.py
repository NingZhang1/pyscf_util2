# coding=UTF-8
import pyscf
import os
import sys
import numpy
import struct
from pyscf import tools
import copy
from pyscf.lib import logger
from pyscf import lib
from pyscf import ao2mo
from pyscf import mcscf, fci, gto, scf

from pyscf_util.MeanField.iciscf import iCI

BASIS = 'ccpvtz'

mol = gto.M(
    verbose=4,
    atom="""
            C   0.000000000000       0.000000000000      -0.621265
            O   0.000000000000       0.000000000000       0.621265
            """,
    basis=BASIS,
    spin=0,
    charge=0,
    symmetry="c2v",
    unit="angstrom"
)
mol.build()
mf = scf.RHF(mol)
mf.kernel()

norb = 9
nelec = 12
mymc2step = mcscf.CASSCF(mf, norb, nelec)
mymc2step.fcisolver = iCI(
    mol=mol,
    cmin=0.0,
    state=[[0, 0, 1]],
    tol=1e-12,
    mo_coeff=mf.mo_coeff,
    taskname="iCI0",
    CVS=True
)
mymc2step.fcisolver.config["segment"] = "0 0 1 4 4 0 0 0"
mymc2step.fcisolver.config["selection"] = 1
mymc2step.fcisolver.config["nvalelec"] = 12
mymc2step.mc1step()

# dump fcidump #

mf.mo_coeff = mymc2step.mo_coeff
tools.fcidump.from_scf(mf, "FCIDUMP_CO_C_Kedge_%s"%BASIS, tol=1e-10)