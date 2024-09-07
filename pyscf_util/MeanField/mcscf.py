import pyscf
from pyscf import fci, mcscf
import numpy as np


def kernel(
    _mol,
    _rohf,
    # _mo_init,
    _nelecas,
    _ncas,
    _frozen=None,
    _mo_init=None,
    _cas_list=None,
    _mc_conv_tol=1e-7,
    _mc_max_macro=128,
    _pyscf_state=None,  # [spintwo, irrep, nstates]
    _do_pyscf_analysis=False,
    _run_mcscf=True,
):
    """Run MCSCF calculation

    TODO: add iCI solver

    Args:
        _mol                : a molecule object
        _rohf               : a pyscf ROHF object
        _nelecas            : a tuple (nelec_alpha, nelec_beta) or an integer (= nelec_alpha + nelec_beta)
        _ncas               : the number of active orbitals
        _frozen             : the number of frozen orbitals
        _mo_init            : the initial guess of MO coefficients
        _cas_list           : the list of active orbitals
        _mc_conv_tol        : the convergence tolerance of MCSCF
        _mc_max_macro       : the maximum number of macro iterations
        _pyscf_state        : the list of states to be calculated
        _do_pyscf_analysis  : whether to do pyscf analysis
        _internal_rotation  : whether to use internal rotation
        _run_mcscf          : whether to run MCSCF

    Kwargs:

    Returns:
        my_mc     : a pyscf MCSCF object
        (mo_init) : the initial guess of MO coefficients

    """

    # Generate MCSCF object
    my_mc = pyscf.mcscf.CASSCF(_rohf, nelecas=_nelecas, ncas=_ncas)
    my_mc.conv_tol = _mc_conv_tol
    my_mc.max_cycle_macro = _mc_max_macro
    my_mc.frozen = _frozen
    # Sort MO
    if _mo_init is None:
        mo_init = _rohf.mo_coeff
    else:
        mo_init = _mo_init
    my_mc.mo_coeff = mo_init  # in case _run_mcscf = False,
    if _cas_list is not None:
        if isinstance(_cas_list, dict):
            mo_init = pyscf.mcscf.sort_mo_by_irrep(my_mc, mo_init, _cas_list)
        else:
            mo_init = pyscf.mcscf.sort_mo(my_mc, mo_init, _cas_list)
    # determine FCIsolver
    nelecas = _nelecas
    if isinstance(_nelecas, tuple):
        nelecas = _nelecas[0] + _nelecas[1]
    if _ncas > 12 or nelecas > 12:  # Large Cas Space
        raise ValueError("Large Cas Space")
    if _pyscf_state is not None:
        # [spintwo, irrep, nstates]
        # Only spin averaged is supported
        solver_all = []
        nstates = 0
        for state in _pyscf_state:
            print(state)
            if state[0] % 2 == 1:
                print("spin 1 is constructed")
                solver = fci.direct_spin1_symm.FCI(_mol)
                solver.wfnsym = state[1]
                solver.nroots = state[2]
                solver.spin = state[0]
                solver_all.append(solver)
            else:
                print("spin 0 is constructed")
                solver = fci.direct_spin0_symm.FCI(_mol)
                solver.wfnsym = state[1]
                solver.nroots = state[2]
                solver.spin = state[0]
                solver_all.append(solver)
            nstates += state[2]
        # exit(1)
        my_mc = mcscf.state_average_mix_(
            my_mc, solver_all, (np.ones(nstates) / nstates)
        )
    # Run
    if _run_mcscf:
        my_mc.kernel(mo_init)
    # Analysis
    if _do_pyscf_analysis:
        my_mc.analyze()
    # 分析成分
    # ao_labels, ao_basis = _get_unique_aolabels(_mol)
    # analysis_mo_contribution(_mol, ao_labels, ao_basis, my_mc.mo_coeff)
    if _run_mcscf:
        return my_mc
    else:
        return [my_mc, mo_init]
