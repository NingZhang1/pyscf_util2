#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PySCF结构优化工具函数
提供便捷的分子结构优化功能
"""

import numpy as np
from pyscf import gto, scf, dft, cc, mcscf

# from pyscf.geomopt import berny_solver, ase_solver
from pyscf.geomopt import berny_solver
from pyscf.geomopt.berny_solver import optimize
from pyscf import lib
import warnings


def geometry_optimization(
    mol_geometry,
    spin=0,
    symmetry=None,
    sfx2c=False,
    charge=0,
    basis="def2-TZVP",
    method="HF",
    max_iter=100,
    conv_tol=1e-8,
    conv_tol_grad=1e-6,
    conv_tol_energy=1e-8,
    solver="berny",
    verbose=4,
    **kwargs,
):
    """
    PySCF分子结构优化函数

    参数:
    ----------
    mol_geometry : str or list
        分子构型，可以是xyz格式的字符串或原子坐标列表
        格式: "原子符号 x y z" 每行一个原子
    spin : int, optional
        自旋多重度 (2S+1)，默认0 (单重态)
    symmetry : str, optional
        分子对称性，如 'C2v', 'D2h' 等，默认None (自动检测)
    sfx2c : bool, optional
        是否使用sf-X2C相对论方法，默认False
    charge : int, optional
        分子电荷，默认0
    basis : str, optional
        基组名称，默认'def2-TZVP'
    method : str, optional
        计算方法，支持: 'HF', 'DFT', 'CCSD', 'CCSD(T)', 'CASSCF'
        默认'HF'
    max_iter : int, optional
        最大优化迭代次数，默认100
    conv_tol : float, optional
        收敛容差，默认1e-6
    conv_tol_grad : float, optional
        梯度收敛容差，默认1e-4
    conv_tol_energy : float, optional
        能量收敛容差，默认1e-6
    solver : str, optional
        优化器类型: 'berny' 或 'ase'，默认'berny'
    verbose : int, optional
        输出详细程度，默认4
    **kwargs : dict
        其他传递给PySCF的参数

    返回:
    ----------
    dict
        包含优化结果的字典:
        - 'mol': 优化后的分子对象
        - 'energy': 最终能量
        - 'coordinates': 优化后的坐标
        - 'converged': 是否收敛
        - 'iterations': 迭代次数
        - 'gradient_norm': 最终梯度范数
    """

    # 设置PySCF参数
    # lib.param.TMPDIR = './tmp'
    # lib.param.VERBOSE = verbose

    # 解析分子构型
    if isinstance(mol_geometry, str):
        # 如果是xyz格式字符串
        mol = parse_xyz_string(mol_geometry)
    elif isinstance(mol_geometry, list):
        # 如果是原子坐标列表
        mol = create_mol_from_coords(mol_geometry)
    else:
        raise ValueError("mol_geometry必须是字符串或列表格式")

    # 设置分子属性
    mol.charge = charge
    mol.spin = spin
    mol.basis = basis
    mol.symmetry = symmetry
    mol.verbose = verbose
    mol.build()

    # 选择计算方法
    if method.upper() == "HF":
        mf = scf.RHF(mol)
    elif method.upper() == "ROHF":
        mf = scf.ROHF(mol)
    elif method.upper() == "DFT":
        functional = kwargs.get("functional", "B3LYP")
        mf = dft.RKS(mol, xc=functional)
    elif method.upper() in ["CCSD", "CCSD(T)"]:
        mf = scf.RHF(mol)
        mf.run()
        if method.upper() == "CCSD":
            mf = cc.CCSD(mf)
        else:
            mf = cc.CCSD(T)(mf)
    elif method.upper() == "CASSCF":
        ncas = kwargs.get("ncas", 2)
        nelecas = kwargs.get("nelecas", 2)
        mf = scf.RHF(mol)
        mf.run()
        mf = mcscf.CASSCF(mf, ncas, nelecas)
    else:
        raise ValueError(f"不支持的方法: {method}")

    # 如果使用sf-X2C
    if sfx2c:
        try:
            from pyscf import x2c

            mf = x2c.sfx2c1e(mf)
        except ImportError:
            warnings.warn("sf-X2C模块未安装，将使用非相对论方法")

    # 运行初始计算
    try:
        mf.run()
    except Exception as e:
        warnings.warn(f"初始计算失败: {e}")
        # 尝试使用更简单的设置
        mf = scf.RHF(mol)
        mf.run()

    # 选择优化器
    if solver.lower() == "berny":
        optimizer = berny_solver.GeometryOptimizer(mf)
    elif solver.lower() == "ase":
        optimizer = ase_solver.GeometryOptimizer(mf)
    else:
        raise ValueError(f"不支持的优化器: {solver}")

    # 设置优化参数
    optimizer.max_iter = max_iter
    optimizer.conv_tol = conv_tol
    optimizer.conv_tol_grad = conv_tol_grad
    optimizer.conv_tol_energy = conv_tol_energy

    # 执行结构优化
    print(f"开始{method}方法的结构优化...")
    print(f"初始能量: {mf.e_tot:.8f} Hartree")

    try:
        mol_eq = optimizer.optimize()

        # rerun #

        # 选择计算方法
        if method.upper() == "HF":
            mf = scf.RHF(mol_eq)
        elif method.upper() == "ROHF":
            mf = scf.ROHF(mol_eq)
        elif method.upper() == "DFT":
            functional = kwargs.get("functional", "B3LYP")
            mf = dft.RKS(mol_eq, xc=functional)
        elif method.upper() in ["CCSD", "CCSD(T)"]:
            mf = scf.RHF(mol_eq)
            mf.run()
            if method.upper() == "CCSD":
                mf = cc.CCSD(mf)
            else:
                mf = cc.CCSD(T)(mf)
        elif method.upper() == "CASSCF":
            ncas = kwargs.get("ncas", 2)
            nelecas = kwargs.get("nelecas", 2)
            mf = scf.RHF(mol_eq)
            mf.run()
            mf = mcscf.CASSCF(mf, ncas, nelecas)
        else:
            raise ValueError(f"不支持的方法: {method}")

        # 如果使用sf-X2C
        if sfx2c:
            try:
                from pyscf import x2c

                mf = x2c.sfx2c1e(mf)
            except ImportError:
                warnings.warn("sf-X2C模块未安装，将使用非相对论方法")

        # rerun #
        mf.run()

        print(f"优化完成!")
        print(f"最终能量: {mf.e_tot:.8f} Hartree")

        # 获取优化结果
        result = {
            "mol": mol_eq,
            "energy": mf.e_tot,
            "coordinates": mol_eq.atom_coords(),
            "converged": optimizer.converged,
            "mf": mf,
            # 'iterations': optimizer.iter_count,
            # 'gradient_norm': optimizer.grad_norm
        }

        return result

    except Exception as e:
        print(f"结构优化失败: {e}")
        return None


def parse_xyz_string(xyz_string):
    """
    解析xyz格式字符串并创建分子对象

    参数:
    ----------
    xyz_string : str
        xyz格式的分子构型字符串

    返回:
    ----------
    pyscf.gto.Mole
        分子对象
    """
    lines = xyz_string.strip().split("\n")

    # 跳过注释行和空行
    atoms = []
    for line in lines:
        line = line.strip()
        if line and not line.startswith("#"):
            parts = line.split()
            if len(parts) >= 4:
                atom_symbol = parts[0]
                x, y, z = map(float, parts[1:4])
                atoms.append([atom_symbol, [x, y, z]])

    # 创建分子对象
    mol = gto.Mole()
    mol.atom = atoms
    return mol


def create_mol_from_coords(atom_list):
    """
    从原子坐标列表创建分子对象

    参数:
    ----------
    atom_list : list
        原子坐标列表，格式: [['原子符号', [x, y, z]], ...]

    返回:
    ----------
    pyscf.gto.Mole
        分子对象
    """
    mol = gto.Mole()
    mol.atom = atom_list
    return mol


def print_optimization_summary(result):
    """
    打印结构优化结果摘要

    参数:
    ----------
    result : dict
        结构优化的结果字典
    """
    if result is None:
        print("结构优化失败")
        return

    print("\n" + "=" * 50)
    print("结构优化结果摘要")
    print("=" * 50)
    print(f"最终能量: {result['energy']:.8f} Hartree")
    print(f"收敛状态: {'是' if result['converged'] else '否'}")
    # print(f"迭代次数: {result['iterations']}")
    # print(f"梯度范数: {result['gradient_norm']:.2e}")

    print("\n优化后的分子构型:")
    mol = result["mol"]
    # for i, (atom, coord) in enumerate(zip(mol.atom_symbol(), mol.atom_coords())):
    #     print(f"{atom:2s} {coord[0]:12.6f} {coord[1]:12.6f} {coord[2]:12.6f}")

    for i, coord in enumerate(mol.atom_coords()):
        atom = mol.atom_symbol(i)
        print(f"{atom:2s} {i+1:2d} {coord[0]:12.8f} {coord[1]:12.8f} {coord[2]:12.8f}")

    # print(mol.atom_symbol())
    # print(mol.atom_coords())
    # exit(1)


# 示例用法
if __name__ == "__main__":
    # 示例1: 水分子优化
    water_xyz = """
    O  0.000000  0.000000  0.000000
    H  0.000000  0.000000  1.000000
    H  0.000000  1.000000  0.000000
    """

    print("示例1: 水分子HF结构优化")
    result1 = geometry_optimization(
        mol_geometry=water_xyz,
        spin=0,
        method="HF",
        basis="cc-pvdz",
        verbose=3,
        symmetry="C2v",
    )
    print_optimization_summary(result1)

    # 示例2: 使用DFT方法
    print("\n示例2: 水分子DFT结构优化")
    result2 = geometry_optimization(
        mol_geometry=water_xyz,
        method="DFT",
        functional="B3LYP",
        basis="6-31G",
        verbose=3,
    )
    print_optimization_summary(result2)
