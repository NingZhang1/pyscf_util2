import numpy
from pyscf.tools import mo_mapping
from pyscf.data.elements import _std_symbol


# def core_orb_id(mol, mo_coeff, atm_label, orb_type="1s"):
#     orb_label = "%s %s" % (atm_label, orb_type)
#     comp = mo_mapping.mo_comps(orb_label, mol, mo_coeff)
#     tmp = list(enumerate(numpy.argsort(-comp)))
#     return tmp[0][1]


# def core_orb_id_dooh(mol, mo_coeff, atm_label, orb_type="1s"):
#     orb_label = "%s %s" % (atm_label, orb_type)
#     comp = mo_mapping.mo_comps(orb_label, mol, mo_coeff)
#     tmp = list(enumerate(numpy.argsort(-comp)))
#     return tmp[0][1], tmp[1][1]


# def sigma_orb_id(mol, mo_coeff, atm_label, orb_type="2pz"):
#     orb_label = "%s %s" % (atm_label, orb_type)
#     comp = mo_mapping.mo_comps(orb_label, mol, mo_coeff)
#     tmp = list(enumerate(numpy.argsort(-comp)))
#     sigma_id, sigma_ast_id = tmp[0][1], tmp[1][1]
#     if sigma_id > sigma_ast_id:
#         sigma_id, sigma_ast_id = sigma_ast_id, sigma_id
#     return sigma_id, sigma_ast_id


# def pi_orb_id(mol, mo_coeff, atm_label, orb_type=["2px", "2py"]):
#     assert mol.natm == 2
#     # for px #
#     orb_label = "%s %s" % (atm_label, orb_type[0])
#     comp = mo_mapping.mo_comps(orb_label, mol, mo_coeff)
#     tmp = list(enumerate(numpy.argsort(-comp)))
#     pi_x_id, pi_ast_x_id = tmp[0][1], tmp[1][1]
#     if pi_x_id > pi_ast_x_id:
#         pi_x_id, pi_ast_x_id = pi_ast_x_id, pi_x_id
#     # for py #
#     orb_label = "%s %s" % (atm_label, orb_type[1])
#     comp = mo_mapping.mo_comps(orb_label, mol, mo_coeff)
#     tmp = list(enumerate(numpy.argsort(-comp)))
#     pi_y_id, pi_ast_y_id = tmp[0][1], tmp[1][1]
#     if pi_y_id > pi_ast_y_id:
#         pi_y_id, pi_ast_y_id = pi_ast_y_id, pi_y_id
#     return [pi_x_id, pi_y_id], [pi_ast_x_id, pi_ast_y_id]
