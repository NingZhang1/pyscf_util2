import numpy


def Dump_Relint_csv(TaskName, relint):
    filename = TaskName + ".csv"
    FILE = open(filename, "w")
    FILE.write("type,i,j,integrals\n")
    for i in range(relint.shape[0]):
        for j in range(relint.shape[1]):
            for k in range(relint.shape[2]):
                if abs(relint[i][j][k]) > 1e-8:
                    FILE.write("%d,%d,%d,%20.12e\n" % (i, j, k, relint[i][j][k]))
    FILE.close()


def Dump_Relint_iCI(filename, relint, nao):
    FILE = open(filename, "w")
    for k in [0, 1, 2, 3]:
        for i in range(0, nao):
            for j in range(0, nao):
                if abs(relint[k][i][j]) > 1e-12:
                    FILE.write(
                        "%.15f %d %d %d\n" % (relint[k][i][j], i, j, k)
                    )  # 0 : x; 1 : y; 2 : z; 3: sf
    FILE.close()


def ReadIn_Relint_csv(TaskName, nao, skiprows=1):
    filename = TaskName + ".csv"
    i, j, k, val = numpy.loadtxt(
        filename,
        dtype=numpy.dtype("i,i,i,d"),
        delimiter=",",
        skiprows=skiprows,
        unpack=True,
    )
    relint = numpy.zeros((4, nao, nao))
    relint[i, j, k] = val
    return relint
