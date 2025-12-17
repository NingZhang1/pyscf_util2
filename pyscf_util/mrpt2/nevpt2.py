#!/usr/bin/env python
#
# Authors: Ning Zhang <ningzhang1024@gmail.com>
#

import ctypes

import tempfile
from functools import reduce
import numpy
import h5py
from pyscf import lib
from pyscf.lib import logger
from pyscf import fci
from pyscf.mcscf import casci, mc1step, mc_ao2mo
from pyscf import ao2mo
from pyscf.ao2mo import _ao2mo
from pyscf import mrpt
import numpy

# pyscf util

from pyscf_util.File import file_rdm

libmc = lib.load_library("libmcscf")

NUMERICAL_ZERO = 1e-14

# 将 iCI rdm 转成 pyscf rdm12

# ################
# Utils
# ################


def _read_fock(fock_filename, norb):
    fock = numpy.zeros((norb, norb))

    i, j, val = numpy.loadtxt(
        fock_filename,
        dtype=numpy.dtype("i,i,d"),
        delimiter=",",
        skiprows=1,
        unpack=True,
    )

    fock[i, j] = fock[j, i] = val

    return fock


def _make_pyscf_rdm12(
    file_rdm1: str, file_rdm2: str, ncas, with_core=False, ncore=0  #
):
    if with_core:
        norb = ncas + ncore
    else:
        norb = ncas

    # ############
    # readin rdm1
    # ############

    i, j, val = numpy.loadtxt(
        file_rdm1, dtype=numpy.dtype("i,i,d"), delimiter=",", skiprows=1, unpack=True
    )

    rdm1 = numpy.zeros((norb, norb))
    rdm1[i, j] = rdm1[j, i] = val

    # remove core orbitals

    if with_core:
        rdm1 = rdm1[ncore:, ncore:]

    # ############
    # readin rdm2
    # ############

    i, j, k, l, val = numpy.loadtxt(
        file_rdm2,
        dtype=numpy.dtype("i,i,i,i,d"),
        delimiter=",",
        skiprows=1,
        unpack=True,
    )

    rdm2 = numpy.zeros((norb, norb, norb, norb))
    rdm2[i, j, k, l] = rdm2[j, i, l, k] = val
    rdm2 = rdm2.transpose(0, 3, 1, 2)  # p^+ q r^+ s

    # remove core orbitals

    if with_core:
        rdm2 = rdm2[ncore:, ncore:, ncore:, ncore:]

    # ####################################
    # in mrpt2 module rdm2 = <p^+q r^+s>
    # we have to add the contr from rdm1
    # ####################################

    for i in range(norb):
        rdm2[:, i, i, :] += rdm1

    return rdm1, rdm2


# make ERIS


def _make_ERIS(mc, mo_coeff, method="incore"):
    from pyscf.mrpt.nevpt2 import _ERIS

    return _ERIS(mc, mo_coeff, method=method)


# ################
# SC-NEVPT2
# ################

from pyscf.mrpt.nevpt2 import Sijrs, Sijr, Srsi


def sc_nevpt2_ici(
    mc,
    eris,
    file_rdm1,
    file_rdm2,
    with_core=False,
    eris_method="incore",
    fock_filename=None,
    cmoao_filename=None,
):

    # ########################
    # readin rdm1 and rdm2
    # ########################

    ncas = mc.ncas
    ncore = mc.ncore

    rdm1, rdm2 = _make_pyscf_rdm12(
        file_rdm1, file_rdm2, ncas, with_core=with_core, ncore=ncore
    )

    dms = {
        "1": rdm1,
        "2": rdm2,
    }

    # ########################
    # readin cmoao
    # ########################

    nmo = mc.mol.nao_nr()

    if cmoao_filename is not None:
        cmoao = numpy.zeros((nmo, nmo))
        i, j, val = numpy.loadtxt(
            cmoao_filename,
            dtype=numpy.dtype("i,i,d"),
            delimiter=",",
            skiprows=1,
            unpack=True,
        )
        cmoao[i, j] = val
        mc.mo_coeff = cmoao

    # ########################

    # ########################
    # make ERIS
    # ########################

    if eris is None:
        eris = _make_ERIS(mc, mc.mo_coeff, method=eris_method)

    # ########################
    # make fock
    # ########################

    if fock_filename is not None:
        fock = _read_fock(fock_filename, mc.mo_coeff.shape[1])
        # set mo_energy #
        mc.mo_energy = numpy.diag(fock)

    # ########################
    # do ijkl
    # ########################

    norm_ijrs, e_ijrs = Sijrs(mc, eris)

    # ########################
    # do ij,a
    # ########################

    norm_ijr, e_ijr = Sijr(mc, dms, eris)

    # ########################
    # do ab,i
    # ########################

    norm_rsi, e_rsi = Srsi(mc, dms, eris)

    res = {
        "ijrs": {
            "norm": norm_ijrs,
            "e": e_ijrs,
        },
        "ijr": {
            "norm": norm_ijr,
            "e": e_ijr,
        },
        "rsi": {
            "norm": norm_rsi,
            "e": e_rsi,
        },
    }

    return res, eris


# ################
# PC-NEVPT2
# ################

from pyscf.mrpt.nevpt2 import make_k27, _extract_orbs, make_hdm1, make_a3
import scipy.linalg
import numpy as np

# def _diag_K_R_rsi(h1e, h2e, dm1, dm2):
#     k27 = make_k27(h1e, h2e, dm1, dm2)
#     eigvals, eigvecs = scipy.linalg.eigh(k27, dm1)
#     return eigvals, eigvecs


def _diag_K_R_rsi(h1e, h2e, dm1, dm2):
    """
    求解广义特征值问题 K * c = e * dm1 * c

    参数:
        h1e, h2e: 哈密顿量矩阵元
        dm1, dm2: 密度矩阵
        tol: 奇异值容忍度

    返回:
        eigvals: 特征值
        eigvecs: 特征向量
    """
    k27 = make_k27(h1e, h2e, dm1, dm2)
    n = dm1.shape[0]
    # 方法1: 使用伪逆处理奇异矩阵
    try:
        # 尝试直接求解广义特征值问题
        eigvals, eigvecs = scipy.linalg.eigh(k27, dm1)
        return eigvals, eigvecs
    except scipy.linalg.LinAlgError:
        # 如果直接求解失败，使用正则化方法
        dm1_regularized = dm1 + 1e-10 * np.eye(n)
        eigvals, eigvecs = scipy.linalg.eigh(k27, dm1_regularized)
        return eigvals, eigvecs


def Srsi_PC(mc, dms, eris, verbose=None):
    # Subspace S_ijr^{(1)}
    mo_core, mo_cas, mo_virt = _extract_orbs(mc, mc.mo_coeff)
    dm1 = dms["1"]
    dm2 = dms["2"]
    ncore = mo_core.shape[1]
    ncas = mo_cas.shape[1]
    nocc = ncore + ncas
    nvirt = mo_virt.shape[1]
    if eris is None:
        h1e = mc.h1e_for_cas()[0]
        h2e = ao2mo.restore(1, mc.ao2mo(mo_cas), ncas).transpose(0, 2, 1, 3)
        h2e_v = ao2mo.incore.general(
            mc._scf._eri, [mo_virt, mo_core, mo_virt, mo_cas], compact=False
        )
        h2e_v = h2e_v.reshape(
            mo_virt.shape[1], ncore, mo_virt.shape[1], ncas
        ).transpose(0, 2, 1, 3)
    else:
        h1e = eris["h1eff"][ncore:nocc, ncore:nocc]
        h2e = eris["ppaa"][ncore:nocc, ncore:nocc].transpose(0, 2, 1, 3)
        h2e_v = eris["pacv"][nocc:].transpose(3, 0, 2, 1)

    e, c = _diag_K_R_rsi(h1e, h2e, dm1, dm2)  # e is the -e)\mu in eqn(31)

    # build S matrix
    S = numpy.einsum("ij,il->lj", c, dm1)
    # make (rs,\mu i) and (rs,\mui')

    tmp1 = numpy.einsum("rsip,pq->rsiq", h2e_v, S)
    rsi_mu_1 = (tmp1 + tmp1.transpose(1, 0, 2, 3)) * (1.0 / numpy.sqrt(2.0))
    rsi_mu_2 = (tmp1 - tmp1.transpose(1, 0, 2, 3)) * (numpy.sqrt(1.5))

    # numerator

    vi_diag = numpy.diag_indices(nvirt)
    vi_triu = numpy.triu_indices(nvirt)

    numerator = -(rsi_mu_1**2 + rsi_mu_2**2)
    # for i in range(numerator.shape[0]):
    numerator[vi_diag] *= 0.5

    denominator = (
        mc.mo_energy[nocc:, None, None, None]
        + mc.mo_energy[None, nocc:, None, None]
        - mc.mo_energy[None, None, :ncore, None]
        + e[None, None, None, :]
    )

    # we only need r<=s

    numerator = numerator[vi_triu]
    denominator = denominator[vi_triu]

    # denominator

    e_rsi = numpy.sum(numerator / denominator)

    return e_rsi

    # print(e)
    # print(c)


def _diag_K_R_ijr(h1e, h2e, dm1, dm2, hdm1):

    a3 = make_a3(h1e, h2e, dm1, dm2, hdm1)
    n = hdm1.shape[0]

    try:
        eigvals, eigvecs = scipy.linalg.eigh(a3, hdm1)
        return eigvals, eigvecs
    except scipy.linalg.LinAlgError:
        # 如果直接求解失败，使用正则化方法
        hdm1_regularized = hdm1 + 1e-10 * np.eye(n)
        eigvals, eigvecs = scipy.linalg.eigh(a3, hdm1_regularized)
        return eigvals, eigvecs


def Sijr_PC(mc, dms, eris, verbose=None):
    # Subspace S_ijr^{(1)}
    mo_core, mo_cas, mo_virt = _extract_orbs(mc, mc.mo_coeff)
    dm1 = dms["1"]
    dm2 = dms["2"]
    ncore = mo_core.shape[1]
    ncas = mo_cas.shape[1]
    nocc = ncore + ncas
    if eris is None:
        h1e = mc.h1e_for_cas()[0]
        h2e = ao2mo.restore(1, mc.ao2mo(mo_cas), ncas).transpose(0, 2, 1, 3)
        h2e_v = ao2mo.incore.general(
            mc._scf._eri, [mo_virt, mo_core, mo_cas, mo_core], compact=False
        )
        h2e_v = h2e_v.reshape(mo_virt.shape[1], ncore, ncas, ncore).transpose(
            0, 2, 1, 3
        )
    else:
        h1e = eris["h1eff"][ncore:nocc, ncore:nocc]
        h2e = eris["ppaa"][ncore:nocc, ncore:nocc].transpose(0, 2, 1, 3)
        h2e_v = eris["pacv"][:ncore].transpose(3, 1, 2, 0)
    if "h1" in dms:
        hdm1 = dms["h1"]
    else:
        hdm1 = make_hdm1(dm1)

    e, c = _diag_K_R_ijr(h1e, h2e, dm1, dm2, hdm1)

    # build S matrix
    S = numpy.einsum("ij,il->lj", c, hdm1)

    # make (r\mu,ji) and (r\mu',ji)
    tmp1 = numpy.einsum("raji,ap->rpji", h2e_v, S)
    r_mu_ji_1 = (tmp1 + tmp1.transpose(0, 1, 3, 2)) * (1.0 / numpy.sqrt(2.0))
    r_mu_ji_2 = (tmp1 - tmp1.transpose(0, 1, 3, 2)) * (numpy.sqrt(1.5))

    # numerator
    numerator = -(r_mu_ji_1**2 + r_mu_ji_2**2)
    # for i in range(numerator.shape[0]):
    ci_diag = numpy.diag_indices(ncore)
    ci_triu = numpy.triu_indices(ncore)

    numerator[:, :, ci_diag[0], ci_diag[1]] *= 0.5

    # denominator
    denominator = (
        mc.mo_energy[nocc:, None, None, None]
        + e[None, :, None, None]
        - mc.mo_energy[None, None, :ncore, None]
        - mc.mo_energy[None, None, None, :ncore]
    )

    numerator = numerator[:, :, ci_triu[0], ci_triu[1]]
    denominator = denominator[:, :, ci_triu[0], ci_triu[1]]

    e_ijr = numpy.sum(numerator / denominator)

    return e_ijr


# ################
# PC driver
# ################


def pc_nevpt2_ici(
    mc,
    eris,
    file_rdm1,
    file_rdm2,
    with_core=False,
    eris_method="incore",
    fock_filename=None,
    cmoao_filename=None,
):

    # ########################
    # readin rdm1 and rdm2
    # ########################

    ncas = mc.ncas
    ncore = mc.ncore

    rdm1, rdm2 = _make_pyscf_rdm12(
        file_rdm1, file_rdm2, ncas, with_core=with_core, ncore=ncore
    )

    dms = {
        "1": rdm1,
        "2": rdm2,
    }

    # ########################
    # readin cmoao
    # ########################

    nmo = mc.mol.nao_nr()

    if cmoao_filename is not None:
        cmoao = numpy.zeros((nmo, nmo))
        i, j, val = numpy.loadtxt(
            cmoao_filename,
            dtype=numpy.dtype("i,i,d"),
            delimiter=",",
            skiprows=1,
            unpack=True,
        )
        cmoao[i, j] = val
        mc.mo_coeff = cmoao

    # ########################
    # make ERIS
    # ########################

    if eris is None:
        eris = _make_ERIS(mc, mc.mo_coeff, method=eris_method)

    # ########################
    # make fock
    # ########################

    if fock_filename is not None:
        fock = _read_fock(fock_filename, mc.mo_coeff.shape[1])
        # set mo_energy #
        mc.mo_energy = numpy.diag(fock)

    # ########################
    # do ijrs
    # ########################

    _, e_ijrs = Sijrs(mc, eris)

    # ########################
    # do ij,a
    # ########################

    e_ijr = Sijr_PC(mc, dms, eris)

    # ########################
    # do ab,i
    # ########################

    e_rsi = Srsi_PC(mc, dms, eris)

    res = {
        "ijrs": {
            "e": e_ijrs,
        },
        "ijr": {
            "e": e_ijr,
        },
        "rsi": {
            "e": e_rsi,
        },
    }

    return res, eris


if __name__ == "__main__":

    from pyscf import gto, scf, mcscf, fci
    from pyscf_util.Integrals import integral_CASCI
    from pyscf_util.iCIPT2.iCIPT2 import kernel

    ## use Cr2 as the example ##

    b = 0.621265 * 2
    mol = gto.Mole()
    mol.build(
        verbose=10,
        output=None,
        atom=[
            ["C", (0.000000, 0.000000, -b / 2)],
            ["C", (0.000000, 0.000000, b / 2)],
        ],
        basis={
            "C": "cc-pvtz",
        },
        symmetry="d2h",
    )

    # Create HF molecule

    mf = scf.sfx2c(scf.RHF(mol)).run()

    norb = 8
    nelec = 8

    mc = mcscf.CASSCF(mf, norb, nelec)
    solver1 = fci.direct_spin1_symm.FCI(mol)
    solver1.wfnsym = "ag"
    solver1.nroots = 1
    solver1.spin = 0

    mc.mc1step()
    # dump cmoao
    # readin cmoao

    mrpt.nevpt2.sc_nevpt(mc)

    #####################
    # iCIPT2-NEVPT2
    #####################

    mo_coeff = mc.mo_coeff
    ncore = mc.ncore

    # dump heff

    integral_CASCI.dump_heff_casci(
        mol,
        mc,
        mo_coeff[:, :ncore],
        mo_coeff[:, ncore : ncore + norb],
        _filename="FCIDUMP_C2",
    )

    kernel(
        IsCSF=True,
        task_name="c2_rdm12",
        fcidump="FCIDUMP_C2",
        segment="0 0 4 4 0 0",
        nelec_val=8,
        rotatemo=0,
        cmin=0.0,
        perturbation=0,
        dumprdm=2,
        relative=0,
        Task="0 0 1 1",
        inputocfg=0,
        etol=1e-10,
        selection=1,
        doublegroup=None,
        direct=None,
        start_with=None,
        end_with=[".csv"],
    )

    # run iCI to dump rdm1 and rdm2

    file_rdm1 = "rdm1.csv"
    file_rdm2 = "rdm2.csv"

    # run iCIPT2-NEVPT2

    res, eris = sc_nevpt2_ici(mc, None, file_rdm1, file_rdm2)

    # print res

    for k, v in res.items():
        print(k)
        print("norm %15.8f e %15.8f" % (v["norm"], v["e"]))
        print("-" * 100)

    # #####################
    # PC-NEVPT2
    # #####################

    # dm1, dm2 = _make_pyscf_rdm12(
    #     file_rdm1, file_rdm2, mc.ncas, with_core=False, ncore=mc.ncore
    # )
    # e_rsi = Srsi_PC(mc, {"1": dm1, "2": dm2}, eris)
    # e_ijr = Sijr_PC(mc, {"1": dm1, "2": dm2}, eris)
    # print(e_rsi)
    # print(e_ijr)

    res, eris = pc_nevpt2_ici(mc, eris, file_rdm1, file_rdm2)

    for k, v in res.items():
        print(k)
        print("e %15.8f" % (v["e"]))
        print("-" * 100)
