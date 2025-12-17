# coding=UTF-8

import pyscf
from pyscf import tools
from pyscf import symm
from pyscf.tools import fcidump
import pyscf.mcscf
from pyscf_util.File import file_rdm
import numpy as np


def OrbSymInfo(Mol, mo_coeff):
    IRREP_MAP = {}
    nsym = len(Mol.irrep_name)
    for i in range(nsym):
        IRREP_MAP[Mol.irrep_name[i]] = i
    # print(IRREP_MAP)

    OrbSym = pyscf.symm.label_orb_symm(Mol, Mol.irrep_name, Mol.symm_orb, mo_coeff)
    IrrepOrb = []
    for i in range(len(OrbSym)):
        IrrepOrb.append(symm.irrep_name2id(Mol.groupname, OrbSym[i]))
    return IrrepOrb


def Canonicalize(Mol, mo_coeff, gfock, ncore, nact, nvir):

    orbsym = OrbSymInfo(Mol, mo_coeff)

    core_sym = orbsym[:ncore]
    act_sym = orbsym[ncore : ncore + nact]
    virtual_sym = orbsym[ncore + nact :]

    # Diagonalize gfock in core, act and virtual spaces #

    core_fock = gfock[:ncore, :ncore]
    act_fock = gfock[ncore : ncore + nact, ncore : ncore + nact]
    virtual_fock = gfock[ncore + nact :, ncore + nact :]

    core_unique_sym = sorted(set(core_sym))
    core_transform = np.zeros((ncore, ncore))

    for sym in core_unique_sym:

        # 找到该对称性的轨道索引
        sym_indices = [i for i, s in enumerate(core_sym) if s == sym]
        block_size = len(sym_indices)

        if block_size > 0:
            # 提取对称性块
            block_fock = core_fock[np.ix_(sym_indices, sym_indices)]

            # 对角化该对称性块
            block_eigvals, block_eigvecs = np.linalg.eigh(block_fock)

            # 将本征向量放入变换矩阵的对应位置
            # for i, idx in enumerate(sym_indices):
            #     core_transform[start_idx:start_idx + block_size, idx] = block_eigvecs[:, i]
            core_transform[np.ix_(sym_indices, sym_indices)] = block_eigvecs

    act_unique_sym = sorted(set(act_sym))
    act_transform = np.zeros((nact, nact))

    for sym in act_unique_sym:

        # 找到该对称性的轨道索引
        sym_indices = [i for i, s in enumerate(act_sym) if s == sym]
        block_size = len(sym_indices)

        if block_size > 0:
            # 提取对称性块
            block_fock = act_fock[np.ix_(sym_indices, sym_indices)]

            # 对角化该对称性块
            block_eigvals, block_eigvecs = np.linalg.eigh(block_fock)

            # 将本征向量放入变换矩阵的对应位置
            act_transform[np.ix_(sym_indices, sym_indices)] = block_eigvecs

    virtual_unique_sym = sorted(set(virtual_sym))
    virtual_transform = np.zeros((nvir, nvir))

    for sym in virtual_unique_sym:
        # 找到该对称性的轨道索引
        sym_indices = [i for i, s in enumerate(virtual_sym) if s == sym]
        block_size = len(sym_indices)

        if block_size > 0:
            # 提取对称性块
            block_fock = virtual_fock[np.ix_(sym_indices, sym_indices)]

            # 对角化该对称性块
            block_eigvals, block_eigvecs = np.linalg.eigh(block_fock)

            # 将本征向量放入变换矩阵的对应位置
            virtual_transform[np.ix_(sym_indices, sym_indices)] = block_eigvecs

    # Update mo_coeff with new eigenvectors

    new_mo_coeff = np.zeros_like(mo_coeff)
    new_mo_coeff[:, :ncore] = mo_coeff[:, :ncore] @ core_transform
    new_mo_coeff[:, ncore : ncore + nact] = (
        mo_coeff[:, ncore : ncore + nact] @ act_transform
    )
    new_mo_coeff[:, ncore + nact :] = mo_coeff[:, ncore + nact :] @ virtual_transform

    # full transform #

    trans_mat_old_2_new = np.eye(ncore + nact + nvir)
    trans_mat_old_2_new[:ncore, :ncore] = core_transform
    trans_mat_old_2_new[ncore : ncore + nact, ncore : ncore + nact] = act_transform
    trans_mat_old_2_new[ncore + nact :, ncore + nact :] = virtual_transform

    # transform gfock #

    new_gfock = trans_mat_old_2_new.T @ gfock @ trans_mat_old_2_new

    orb_ene = np.diag(new_gfock)

    return new_mo_coeff, new_gfock, orb_ene
