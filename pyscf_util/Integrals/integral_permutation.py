import numpy as np
from pyscf_util.misc.misc import _combine4


def kernel(h1e, h2e, norb, map_old_2_new, symmetry=4):

    assert symmetry == 4, "Only symmetry 4 is supported"

    h1e_new = np.zeros((norb, norb))
    h2e_new = np.zeros(len(h2e))

    # loop 1e

    for p in range(norb):
        for q in range(norb):
            h1e_new[map_old_2_new[p], map_old_2_new[q]] = h1e[p, q]

    # loop 2e

    indx = 0

    for p in range(norb):
        for q in range(p + 1):
            for r in range(p + 1):
                end_s = q + 1 if p == r else r + 1
                for s in range(end_s):
                    indx_new = _combine4(*(map_old_2_new[x] for x in (p, q, r, s)))
                    h2e_new[indx_new] = h2e[indx]
                    indx += 1

    return h1e_new, h2e_new
