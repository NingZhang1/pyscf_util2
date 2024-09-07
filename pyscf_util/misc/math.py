import numpy as np
from scipy import linalg


def _orthogonalize(vec, ovlp=None):
    if ovlp is None:
        ovlp = np.eye(vec.shape[0])

    # Solve the generalized eigenvalue problem
    w, v = linalg.eigh(vec.T @ ovlp @ vec)

    # Construct orthogonalized vectors
    return vec @ v / np.sqrt(w)
