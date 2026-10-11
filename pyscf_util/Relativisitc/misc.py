"""Wigner rotation D-matrix for J=1/2"""

import numpy
from scipy.special import factorial


def dmatrix_half(beta):
    """Wigner small-d matrix (in z-y-z convention) for J=1/2"""
    c = numpy.cos(beta / 2)
    s = numpy.sin(beta / 2)
    mat = numpy.array(((c, -s), (s, c)), dtype=numpy.complex128)
    return mat


def Dmatrix_half(alpha, beta, gamma):
    """Wigner rotation D-matrix

    D_{mm'} = <lm|R(alpha,beta,gamma)|lm'>
    alpha, beta, gamma are Euler angles (in z-y-z convention)
    l = 1/2

    Kwargs:
        reorder_p (bool): Whether to put the p functions in the (x,y,z) order.
    """
    d = dmatrix_half(beta)
    ms = numpy.asarray([0.5, -0.5])  # |alpha>, |beta>
    D = numpy.einsum(
        "i,ij,j->ij", numpy.exp(-1j * alpha * ms), d, numpy.exp(-1j * gamma * ms)
    )

    return D


if __name__ == "__main__":

    from pyscf.symm.Dmatrix import Dmatrix

    c1 = [1, 0, 0]
    c2 = [-1, 0, 0]

    Ebar = [0, 0, 2 * numpy.pi]
    RX = [numpy.pi, numpy.pi, 0]
    RY = [0, numpy.pi, 0]
    RZ = [0, 0, numpy.pi]

    print(Dmatrix_half(Ebar[0], Ebar[1], Ebar[2]))
    print(Dmatrix_half(RX[0], RX[1], RX[2]))
    print(Dmatrix_half(RY[0], RY[1], RY[2]))
    print(Dmatrix_half(RZ[0], RZ[1], RZ[2]))

    print(Dmatrix(1, RX[0], RX[1], RX[2], True))
    print(Dmatrix(1, RY[0], RY[1], RY[2], True))
    print(Dmatrix(1, RZ[0], RZ[1], RZ[2], True))

    print(numpy.diag(Dmatrix(2, RX[0], RX[1], RX[2], True)))
    print(numpy.diag(Dmatrix(2, RY[0], RY[1], RY[2], True)))
    print(numpy.diag(Dmatrix(2, RZ[0], RZ[1], RZ[2], True)))

    print(numpy.diag(Dmatrix(3, RX[0], RX[1], RX[2], True)))
    print(numpy.diag(Dmatrix(3, RY[0], RY[1], RY[2], True)))
    print(numpy.diag(Dmatrix(3, RZ[0], RZ[1], RZ[2], True)))

    print(numpy.diag(Dmatrix(4, RX[0], RX[1], RX[2], True)))
    print(numpy.diag(Dmatrix(4, RY[0], RY[1], RY[2], True)))
    print(numpy.diag(Dmatrix(4, RZ[0], RZ[1], RZ[2], True)))
