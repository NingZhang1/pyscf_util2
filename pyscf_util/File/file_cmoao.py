import numpy


def Dump_Cmoao(TaskName, mocoeff):
    filename = TaskName + ".csv"
    FILE = open(filename, "w")
    FILE.write("i,j,mocoeff\n")
    for i in range(mocoeff.shape[0]):
        for j in range(mocoeff.shape[1]):
            FILE.write("%d,%d,%20.12e\n" % (i, j, mocoeff[i][j]))
    FILE.close()


def ReadIn_Cmoao(TaskName, nao, nmo=None, skiprows=1):
    filename = TaskName + ".csv"
    i, j, val = numpy.loadtxt(
        filename,
        dtype=numpy.dtype("i,i,d"),
        delimiter=",",
        skiprows=skiprows,
        unpack=True,
    )
    if nmo is None:
        nmo = nao
    cmoao = numpy.zeros((nao, nmo))
    cmoao[i, j] = val
    return cmoao
