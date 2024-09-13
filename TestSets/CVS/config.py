GEOMETRY = {
    "CO": """
C 0.000000 0.000000  0.000000
O 0.000000 0.000000  1.128400
""",
    "CO2": """
C 0.000000 0.000000  0.000000
O 0.000000 0.000000  1.160000
O 0.000000 0.000000 -1.160000
""",
    "CH4": """
C   0.0000000  0.0000000  0.0000000
H   0.0000000  0.0000000  1.0861000
H   0.0000000 -1.0239849 -0.3620333
H  -0.8867969  0.5119924 -0.3620333
H   0.8867969  0.5119924 -0.3620333
""",
    "NH3": """
N   0.00000000 0.46869816  0.37922759
H   0.81180903 0.00000000  0.00000000
H   0.00000000 1.40609449  0.00000000
H  -0.81180903 0.00000000  0.00000000
""",
    "H2O": """
O  0.000000  0.0000000 -0.1172519
H  0.000000 -0.7571643  0.4690074
H  0.000000  0.7571643  0.4690074
""",
    "HF": """
F  0.000000 0.000000 0.000000
H  0.000000 0.000000 0.916800
""",
    "C2H2": """
C 0.000000 0.000000 -0.601500
C 0.000000 0.000000  0.601500
H 0.000000 0.000000 -1.663000
H 0.000000 0.000000  1.663000
""",
    "H2CO": """
C  0.0000000   0.0000000	 -0.5296279
O  0.0000000   0.0000000	  0.6741721
H  0.0000000   0.9361475	 -1.1078044
H  0.0000000  -0.9361475     -1.1078044
""",
    "N2": """
N 0.000000 0.000000 0.000000
N 0.000000 0.000000 1.097600
""",
    "O2": """
O 0.000000 0.000000 0.000000
O 0.000000 0.000000 1.207700
""",  # NOTE: 基态是个三重态！
    "SiH4": """
Si  0.000000  0.000000  0.000000
H   0.854400  0.854400  0.854400
H  -0.854400 -0.854400  0.854400
H  -0.854400  0.854400 -0.854400
H   0.854400 -0.854400 -0.854400
""",
    "PH3": """
P   0.00000000 0.59312580  0.76512855
H   1.02732402 0.00000000  0.00000000
H   0.00000000 1.77937741  0.00000000
H  -1.02732402 0.00000000  0.00000000
""",
    "H2S": """
S  0.000000  0.000000  -0.1855710
H  0.000000  0.9608222  0.7422841
H  0.000000 -0.9608222  0.7422841
""",
    "HCl": """
Cl  0.000000 0.000000 0.000000
H   0.000000 0.000000 1.274400
""",
    "SO2": """
S  0.0000000    0.0000000    -0.3618316
O  0.0000000    1.2359241     0.3618316
O  0.0000000   -1.2359241     0.3618316
""",
    "H3COH": """
C  -0.0503000   0.6685000    0.0000000
O  -0.0503000  -0.7585000    0.0000000
H  -1.0807000   1.0417000    0.0000000
H   0.4650000   1.0417000    0.8924000
H   0.4650000   1.0417000   -0.8924000
H   0.8544000  -1.0417000    0.0000000
""",
    "Cl2": """
Cl 0.000000 0.000000 0.000000
Cl 0.000000 0.000000 1.987900
""",
    "NNO": """
    N 0.0000    0.0000    0.0000
    N 1.0256    0.0000    0.0000
    O 2.0000    0.0000    0.0000
""",
    "CH3CN": """
   N  1.2608    0.0000    0.0000 
   C -1.3613    0.0000    0.0000 
   C  0.1006    0.0000    0.0000 
   H -1.7500   -0.8301    0.5974 
   H -1.7501   -0.1022   -1.0175 
   H -1.7500    0.9324    0.4202 
   """,
    "HCN": """
  N -0.5800    0.0000    0.0000 
  C  0.5800    0.0000    0.0000 
  H  1.6450    0.0000    0.0000 
""",
    "O3": """
  O -0.0950   -0.4943    0.0000 
  O  1.1489    0.2152    0.0000 
  O -1.0540    0.2791    0.0000 
""",
}  # NOTE: in Bohr


##### GEOMETRY OPTIMIZED VIA B3LYP aug-ccpvtz in Bohr #####

GEOMETRY_OPTIMIZED = {
    "CO": """
C 0.         0.         0.002402  
O 0.         0.         2.12996496
""",
    "CO2": """
C 0.          0.          0.        
O1 0.          0.          2.19297627
O2 0.          0.         -2.19297627
""",
    "CH4": """
C   4.68249694e-14 -1.01868188e-06 -2.69392358e-07
H1   7.35241212e-14 -1.90735245e-07  2.05654391e+00
H2  -1.30449729e-13 -1.93892934e+00 -6.85514063e-01
H3  -1.67916106e+00  9.69465181e-01 -6.85514693e-01
H4   1.67916106e+00  9.69465181e-01 -6.85514693e-01
""",
    "NH3": """
N  -1.45354658e-14  8.85711161e-01  7.10648982e-01
H1   1.54081995e+00 -3.88178152e-03  1.99576692e-03
H2   3.56290566e-14  2.66489705e+00  1.99576837e-03
H3  -1.54081995e+00 -3.88178152e-03  1.99576692e-03
""",
    "H2O": """
O  0.          0.         -0.22000133
H1  0.         -1.44281777  0.88550921
H2  0.          1.44281777  0.88550921
""",
    "HF": """
F  0.          0.         -0.00692527
H  0.          0.          1.73942619
""",
    "C2H2": """
C1 0.          0.         -1.13045612
C2 0.          0.          1.13045612
H1 0.          0.         -3.13661966
H2 0.          0.          3.13661966
""",
    "H2CO": """
C  0.          0.         -0.99260976
O  0.          0.          1.27567854
H1  0.          1.77341306 -2.09840683
H2  0.         -1.77341306 -2.09840683
""",
    "N2": """
N1 0.         0.        -1.0309951
N2 0.         0.         1.0309951
""",
    "O2": """
O1 0.          0.         -1.13791786
O2 0.          0.          1.13791786
""",  # NOTE: 基态是个三重态！
    "SiH4": """
Si -1.53358317e-13  4.83937748e-13 -4.36697869e-13
H1   1.61747223e+00  1.61747223e+00  1.61747223e+00
H2  -1.61747223e+00 -1.61747223e+00  1.61747223e+00
H3  -1.61747223e+00  1.61747223e+00 -1.61747223e+00
H4   1.61747223e+00 -1.61747223e+00 -1.61747223e+00
""",
    "PH3": """
P   8.31537687e-14  1.12084654e+00  1.44917485e+00
H1   1.95612313e+00 -8.52643878e-03 -1.09740995e-03
H2  -1.53388980e-13  3.37958764e+00 -1.09662522e-03
H3  -1.95612313e+00 -8.52643878e-03 -1.09740995e-03
""",
    "H2S": """
S  0.          0.         -0.35232488
H1  0.          1.83662056  1.40353691
H2  0.         -1.83662056  1.40353691
""",
    "HCl": """
Cl  0.          0.         -0.00868048
H   0.          0.          2.41694745
""",
    "SO2": """
S  0.          0.         -0.71086447
O1  0.          2.35338623  0.69731355
O2  0.         -2.35338623  0.69731355
""",
    "H3COH": """
C  -0.06153042  1.27161433  0.        
O  -0.11262059 -1.41864888  0.        
H1  -2.01894885  1.90072611  0.        
H2   0.87002457  2.02742646  1.6835658 
H3   0.87002457  2.02742646 -1.6835658 
H4   1.59274454 -2.04156444  0.        
""",
    "Cl2": """
Cl1 0.          0.         -1.91175365
Cl2 0.          0.          1.91175365
""",
    "NNO": """
    N1 -2.52146955e-01  0.00000000e+00 0.00000000e+00
    N2  1.86622219e+00  0.00000000e+00 0.00000000e+00
    O  4.10348013e+00  0.00000000e+00 0.00000000e+00
""",
    "CH3CN": """
   N  2.34428892e+00 -3.52234270e-06 -2.96776429e-05 
   C1 -2.57636409e+00  3.46414736e-05  4.38525144e-05 
   C2  1.72212306e-01  9.57945287e-06 -2.08898454e-05 
   H1 -3.28705613e+00 -1.56847662e+00  1.12879623e+00 
   H2 -3.28712060e+00 -1.93192614e-01 -1.92267752e+00 
   H3 -3.28702256e+00  1.76181751e+00  7.94076978e-01 
   """,
    "HCN": """
  N -1.07874586e+00  0.00000000e+00 -2.39529697e-16 
  C  1.08686968e+00  0.00000000e+00  2.41333548e-16 
  H  3.10047565e+00  0.00000000e+00  6.88443891e-16 
""",
    "O3": """
  O1 -0.03647225 -0.81030057  0.        
  O2  2.05158405  0.31377831  0.        
  O3 -2.01530078  0.49652226  0.        
""",
}

##############################################################################

# STATE INFO #

# CORE CONFIGURATION #

CORE_CONFIG = {
    "C": ["1s"],
    "N": ["1s"],
    "O": ["1s"],
    "F": ["1s"],
    "Si": ["1s", "2s", "2px", "2py", "2pz"],
    "P": ["1s", "2s", "2px", "2py", "2pz"],
    "S": ["1s", "2s", "2px", "2py", "2pz"],
    "Cl": ["1s", "2s", "2px", "2py", "2pz"],
}

# SIGMA ORB CONFIGURATION #

SIGMA_CONFIG = {
    "C": ["2pz"],
    "N": ["2pz"],
    "O": ["2pz"],
    "F": ["2pz"],
    "Si": ["3pz"],
    "P": ["3pz"],
    "S": ["3pz"],
    "Cl": ["3pz"],
}

PI_CONFIG = {
    "C": ["2px", "2py"],
    "N": ["2px", "2py"],
    "O": ["2px", "2py"],
    "F": ["2px", "2py"],
    "Si": ["3px", "3py"],
    "P": ["3px", "3py"],
    "S": ["3px", "3py"],
    "Cl": ["3px", "3py"],
}

NORB_ANALYSIS = {
    "H": 1,
    "C": 9,
    "N": 9,
    "O": 9,
    "F": 9,
    "Si": 13,
    "P": 13,
    "S": 13,
    "Cl": 13,
}

##############################################################################

# Mole Info #

MOLE_ORB_TYPE = {
    # ---------------------------------------------------------------- #
    "CO": {
        "orb_type": {
            "K_edge_C": [1],
            "K_edge_O": [0],
            "valence": [2, 3, 4, 5, 6, 7, 8, 9],
        },
        "task_type": {
            "K_edge_C": {
                "ordering": ["K_edge_O", "K_edge_C", "valence"],
                "taskinfo": {
                    "all": {"segment": "0 1 1 4 4 0 %d 0"},
                    "fzc": {"segment": "1 0 1 4 4 0 %d 0"},
                    "nleft": 10,
                    "nelec_val": 12,
                },
            },
            "K_edge_O": {
                "ordering": ["K_edge_C", "K_edge_O", "valence"],
                "taskinfo": {
                    "all": {"segment": "0 1 1 4 4 0 %d 0"},
                    "fzc": {"segment": "1 0 1 4 4 0 %d 0"},
                    "nleft": 10,
                    "nelec_val": 12,
                },
            },
        },
    }
    # ---------------------------------------------------------------- #
}

TASK_CASE = ["K_edge", "L_edge", "2s_edge"]
