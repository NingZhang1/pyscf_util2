import numpy


def Dump_SpinRDM1(filename, rdm1):
    norb = rdm1.shape[1]
    nstates = rdm1.shape[0]
    file = open(filename, "w")
    file.write("istate,i,j,rdm1\n")
    for istate in range(nstates):
        for i in range(norb):
            for j in range(norb):
                if abs(rdm1[istate][i][j]) > 1e-10:
                    file.write(
                        "%d,%d,%d,%20.12e\n" % (istate, i, j, rdm1[istate][i][j])
                    )
    file.close()


def Dump_SpinRDM2(filename, rdm2):
    norb = rdm2.shape[1]
    nstates = rdm2.shape[0]
    file = open(filename, "w")
    file.write("istate,p,q,r,s,rdm2\n")
    for istate in range(nstates):
        for p in range(norb):
            for q in range(norb):
                for r in range(norb):
                    for s in range(norb):
                        if abs(rdm2[istate][p][q][r][s]) > 1e-10:
                            file.write(
                                "%d,%d,%d,%d,%d,%20.12e\n"
                                % (istate, p, q, r, s, rdm2[istate][p][q][r][s])
                            )
    file.close()


def ReadIn_SpinRDM1(filename, norb, nstates, IsAveraged=False):
    if IsAveraged:
        i, j, val = numpy.loadtxt(
            filename, dtype=numpy.dtype("i,i,d"), delimiter=",", skiprows=1, unpack=True
        )
        rdm1 = numpy.zeros((norb, norb))
        rdm1[i, j] = val
        return rdm1
    else:
        istate, i, j, val = numpy.loadtxt(
            filename,
            dtype=numpy.dtype("i,i,i,d"),
            delimiter=",",
            skiprows=1,
            unpack=True,
        )
        rdm1 = numpy.zeros((nstates, norb, norb))
        rdm1[istate, i, j] = val
        return rdm1
