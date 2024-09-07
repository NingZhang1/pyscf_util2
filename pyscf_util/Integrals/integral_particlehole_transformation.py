import numpy as np

from pyscf_util.misc.misc import _combine4


def kernel(h1e, h2e, energy_core, norb):

    h1e_new = np.zeros((norb, norb))
    energy_core_new = energy_core

    for i in range(norb):
        for j in range(norb):
            if i == j:
                h1e_new[i, i] = -h1e[i, i] - h2e[_combine4(i, i, i, i)]
                for k in range(norb):
                    if k == i:
                        continue
                    h1e_new[i, i] += (
                        -2 * h2e[_combine4(i, i, k, k)] + h2e[_combine4(i, k, k, i)]
                    )
            else:
                h1e_new[j, i] = (
                    -h1e[i, j] - h2e[_combine4(i, i, i, j)] - h2e[_combine4(i, j, j, j)]
                )
                for k in range(norb):
                    if k == i or k == j:
                        continue
                    h1e_new[j, i] += (
                        h2e[_combine4(i, k, k, j)] - 2 * h2e[_combine4(i, j, k, k)]
                    )

    for i in range(norb):
        energy_core_new += 2 * h1e[i, i] + h2e[_combine4(i, i, i, i)]

    for i in range(norb):
        for j in range(i + 1, norb):
            energy_core_new += (
                4 * h2e[_combine4(i, i, j, j)] - 2 * h2e[_combine4(i, j, j, i)]
            )

    return h1e_new, energy_core_new
