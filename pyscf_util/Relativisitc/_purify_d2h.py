from pyscf import gto, scf, lib, symm
import numpy
from functools import reduce
from pyscf import lib
from pyscf import tools
import pyscf


def label_orb_symm(mol, irrep_name, symm_orb, mo, s=None):
    nmo = mo.shape[1]
    if s is None:
        s = mol.intor_symmetric("int1e_ovlp")
    s_mo = numpy.dot(s, mo)
    norm = numpy.zeros((len(irrep_name), nmo))
    for i, csym in enumerate(symm_orb):
        moso = numpy.dot(csym.T.conj(), s_mo)
        ovlpso = reduce(numpy.dot, (csym.T.conj(), s, csym))
        try:
            s_moso = lib.cho_solve(ovlpso, moso)
        except numpy.linalg.LinAlgError:
            print("numpy.linalg.LinAlgError")
            ovlpso[numpy.diag_indices(csym.shape[1])] += 1e-12
            s_moso = lib.cho_solve(ovlpso, moso)
        norm[i] = numpy.einsum("ki,ki->i", moso.conj(), s_moso).real
    # print(norm)
    norm /= numpy.sum(norm, axis=0)  # for orbitals which are not normalized
    # iridx = numpy.argmax(norm, axis=0)
    # orbsym = numpy.asarray([irrep_name[i] for i in iridx])
    return norm


def purify_pes_d2h(mol, mf):

    # generate ao2L component #

    ao2Lalpha = get_ao2L_spinor(mol)
    ao2Lbeta = get_ao2L_spinor(mol, False)

    # generate mo_Lalpha, mo_Lbeta #

    n2c = mol.nao_2c()
    mo_Lalpha_coeff = ao2Lalpha @ mf.mo_coeff[:n2c, n2c:]
    mo_Lbeta_coeff = ao2Lbeta @ mf.mo_coeff[:n2c, n2c:]
    mo_Lalpha_coeff_real = mo_Lalpha_coeff.real
    mo_Lalpha_coeff_imag = mo_Lalpha_coeff.imag
    mo_Lbeta_coeff_real = mo_Lbeta_coeff.real
    mo_Lbeta_coeff_imag = mo_Lbeta_coeff.imag

    res = mf.mo_coeff[:, n2c:].copy()

    ovlp = mol.intor("int1e_ovlp")

    for i in range(0, n2c):
        norm_alpha_real = (
            mo_Lalpha_coeff_real[:, i].T @ ovlp @ mo_Lalpha_coeff_real[:, i]
        )
        norm_alpha_imag = (
            mo_Lalpha_coeff_imag[:, i].T @ ovlp @ mo_Lalpha_coeff_imag[:, i]
        )
        norm_beta_real = mo_Lbeta_coeff_real[:, i].T @ ovlp @ mo_Lbeta_coeff_real[:, i]
        norm_beta_imag = mo_Lbeta_coeff_imag[:, i].T @ ovlp @ mo_Lbeta_coeff_imag[:, i]

        norm_lst = [norm_alpha_real, norm_alpha_imag, norm_beta_real, norm_beta_imag]
        mo_coeff_lst = numpy.array(
            [
                mo_Lalpha_coeff_real[:, i],
                mo_Lalpha_coeff_imag[:, i],
                mo_Lbeta_coeff_real[:, i],
                mo_Lbeta_coeff_imag[:, i],
            ]
        )

        argmax_loc = numpy.argmax(numpy.array(norm_lst))

        # print(norm_lst)
        # print(argmax_loc)
        # assert argmax_loc in [0, 2]

        symm_loc = label_orb_symm(
            mol,
            mol.irrep_name,
            mol.symm_orb,
            mo_coeff_lst[argmax_loc].reshape(-1, 1),
            ovlp,
            # tol=1e-7/norm_lst[argmax_loc],
        )

        print(norm_lst[argmax_loc])
        print(symm_loc)
        # print(argmax_loc, mol.irrep_name[argmax_loc], norm_lst[argmax_loc])

        # print(
        #     f"MO {i} is assigned to irrep {mol.irrep_name[symm_loc.argmax()]} with norm {norm_lst[argmax_loc]}"
        # )

    # final check #


def purify_pes_c2v(mol, mf):

    # generate ao2L component #

    # #

    return


# purify for d2h and its subgroup #

if __name__ == "__main__":

    from Lalpha_ao2spinor import get_ao2L_spinor

    mol = gto.M(
        atom="""
C   0.000000   0.000000   0.669500
C   0.000000   0.000000  -0.669500
H   0.928900   0.000000   1.232100
H  -0.928900   0.000000   1.232100
H   0.928900   0.000000  -1.232100
H  -0.928900   0.000000  -1.232100
""",
        # basis="unc-cc-pvdz-dk",
        # basis="3-21g",
        basis="cc-pvdz-dk",
        verbose=5,
        charge=0,
        spin=0,
        symmetry="d2h",
    )

    print(get_ao2L_spinor(mol))

    ao2Lalpha = get_ao2L_spinor(mol)
    ao2Lbeta = get_ao2L_spinor(mol, False)

    # do scf #

    mol.build()
    mf = scf.dhf.RDHF(mol)
    mf.conv_tol = 1e-12
    # mf.kernel()
    # purify_pes_d2h(mol, mf)

    # learn symm module line by line #

    sym_org = mol._symm_orig
    sym_axes = mol._symm_axes
    group = mol.groupname

    print("sym_org", sym_org)
    print("sym_axes", sym_axes)
    print("group", group)

    from pyscf.symm import geom, param, basis

    prop_atoms = mol.format_atom(mol._atom, sym_org, sym_axes, "Bohr")
    eql_atom_ids = geom.symm_identical_atoms(group, prop_atoms)

    print(prop_atoms)  # 原子 + 坐标
    print(eql_atom_ids)

    ops = numpy.asarray([param.D2H_OPS[op] for op in param.OPERATOR_TABLE[group]])
    chartab = numpy.array([x[1:] for x in param.CHARACTER_TABLE[group]])
    nirrep = chartab.__len__()
    aoslice = mol.aoslice_by_atom()
    nao = mol.nao_nr()
    atom_coords = numpy.array([a[1] for a in prop_atoms])

    print(ops)  # 以 (x,y,z) 为 basis 的操作矩阵
    print(chartab)
    print(aoslice)
    print(atom_coords)

    print("------------------------------------------------")

    sodic = [[] for i in range(8)]
    for atom_ids in eql_atom_ids:
        print("atom_ids = ", atom_ids)
        r0 = atom_coords[atom_ids[0]]
        print("r0 = ", r0)
        op_coords = numpy.einsum(
            "x,nxy->ny", r0, ops
        )  # ops  作用到第一个原子得到的其他坐标
        print(op_coords)
        # Using ops to generate other atoms from atom_ids[0]
        coords0 = atom_coords[atom_ids]
        print("coords0 = ", coords0)
        print(op_coords.shape, coords0.shape)
        print((op_coords.reshape(-1, 1, 3) - coords0).shape)
        dc = abs(op_coords.reshape(-1, 1, 3) - coords0).max(axis=2)
        print("dc = ", dc)  # 记录操作关联了哪些原子
        op_relate_idx = numpy.argwhere(dc < geom.TOLERANCE)[:, 1]
        print("op_relate_idx = ", op_relate_idx)
        ao_loc = numpy.array([aoslice[atom_ids[i], 2] for i in op_relate_idx])
        print("ao_loc = ", ao_loc)  # 对称操作关联的atom 对应的basis开始的指标

        b0, b1 = aoslice[atom_ids[0], :2]
        ip = 0
        for ib in range(b0, b1):
            l = mol.bas_angular(ib)
            print("ib = ", ib, "l = ", l)
            if mol.cart:  # for spinor 报错 !
                degen = (l + 1) * (l + 2) // 2
                cbase = numpy.zeros((degen, nirrep, nao))
                for op_id, op in enumerate(ops):
                    n = 0
                    for x in range(l, -1, -1):
                        for y in range(l - x, -1, -1):
                            z = l - x - y
                            idx = ao_loc[op_id] + n
                            sign = op[0, 0] ** x * op[1, 1] ** y * op[2, 2] ** z
                            cbase[n, :, idx] += sign * chartab[:, op_id]
                            n += 1
            else:
                degen = l * 2 + 1
                cbase = numpy.zeros((degen, nirrep, nao))
                for op_id, op in enumerate(param.OPERATOR_TABLE[group]):
                    print("op_id = ", op_id, "op = ", op)
                    for n, m in enumerate(range(-l, l + 1)):
                        idx = ao_loc[op_id] + n
                        print("n = ", n, "m = ", m, "idx = ", idx)
                        if basis.tot_parity_odd(op, l, m):
                            cbase[n, :, idx] -= chartab[:, op_id]
                        else:
                            cbase[n, :, idx] += chartab[:, op_id]
            norms = numpy.sqrt(numpy.einsum("mij,mij->mi", cbase, cbase))

            for i in range(mol.bas_nctr(ib)):
                for n, ir in numpy.argwhere(norms > 1e-12):
                    c = numpy.zeros(nao)
                    c[ip:] = cbase[n, ir, : nao - ip] / norms[n, ir]
                    sodic[ir].append(c)
                ip += degen

        print("------------------------------------------------")
