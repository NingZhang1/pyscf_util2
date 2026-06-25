from pyscf import lib
import numpy


from pyscf import __config__

IOBLK_SIZE = getattr(__config__, "ao2mo_outcore_ioblk_size", 256)  # 256 MB
IOBUF_WORDS = getattr(__config__, "ao2mo_outcore_iobuf_words", 1e8)  # 1.6 GB
IOBUF_ROW_MIN = getattr(__config__, "ao2mo_outcore_row_min", 160)
MAX_MEMORY = getattr(__config__, "ao2mo_outcore_max_memory", 4000)  # 4GB

### filename ###


coulomb_LLLL = "%s_coulomb_LLLL.h5"
coulomb_SSSS = "%s_coulomb_SSSS.h5"
coulomb_LLSS = "%s_coulomb_LLSS.h5"
coulomb_SSLL = "%s_coulomb_SSLL.h5"

breit_LSLS = "%s_breit_LSLS.h5"
breit_SLSL = "%s_breit_SLSL.h5"
breit_LSSL = "%s_breit_LSSL.h5"
breit_SLLS = "%s_breit_SLLS.h5"


def _r_outcore_Coulomb(
    mol, my_RDHF, npes, prefix, max_memory=MAX_MEMORY, ioblk_size=IOBLK_SIZE
):

    from pyscf.ao2mo import r_outcore

    n2c = mol.nao_2c()
    mo_coeff = my_RDHF.mo_coeff
    mo_coeff_mat = numpy.matrix(mo_coeff)

    mo_coeff_pes = mo_coeff_mat[:, n2c : n2c + npes]
    mo_coeff_L = mo_coeff_pes[:n2c, :]
    mo_coeff_S = mo_coeff_pes[n2c:, :]

    r_outcore.general(
        mol,
        (mo_coeff_L, mo_coeff_L, mo_coeff_L, mo_coeff_L),
        coulomb_LLLL % prefix,
        intor="int2e_spinor",
        max_memory=max_memory,
        ioblk_size=ioblk_size,
        aosym="s1",
    )
    r_outcore.general(
        mol,
        (mo_coeff_S, mo_coeff_S, mo_coeff_S, mo_coeff_S),
        coulomb_SSSS % prefix,
        intor="int2e_spsp1spsp2_spinor",
        max_memory=max_memory,
        ioblk_size=ioblk_size,
        aosym="s1",
    )
    r_outcore.general(
        mol,
        (mo_coeff_L, mo_coeff_L, mo_coeff_S, mo_coeff_S),
        coulomb_LLSS % prefix,
        intor="int2e_spsp2_spinor",
        max_memory=max_memory,
        ioblk_size=ioblk_size,
        aosym="s1",
    )
    r_outcore.general(
        mol,
        (mo_coeff_S, mo_coeff_S, mo_coeff_L, mo_coeff_L),
        coulomb_SSLL % prefix,
        intor="int2e_spsp1_spinor",
        max_memory=max_memory,
        ioblk_size=ioblk_size,
        aosym="s1",
    )


def _r_outcore_Breit(
    mol, my_RDHF, npes, prefix, max_memory=MAX_MEMORY, ioblk_size=IOBLK_SIZE
):

    from pyscf.ao2mo import r_outcore

    n2c = mol.nao_2c()
    mo_coeff = my_RDHF.mo_coeff
    mo_coeff_mat = numpy.matrix(mo_coeff)

    mo_coeff_pes = mo_coeff_mat[:, n2c : n2c + npes]
    mo_coeff_L = mo_coeff_pes[:n2c, :]
    mo_coeff_S = mo_coeff_pes[n2c:, :]

    r_outcore.general(
        mol,
        (mo_coeff_L, mo_coeff_S, mo_coeff_L, mo_coeff_S),
        breit_LSLS % prefix,
        intor="int2e_breit_ssp1ssp2_spinor",
        max_memory=max_memory,
        ioblk_size=ioblk_size,
        aosym="s1",
    )
    r_outcore.general(
        mol,
        (mo_coeff_S, mo_coeff_L, mo_coeff_S, mo_coeff_L),
        breit_SLSL % prefix,
        intor="int2e_breit_sps1sps2_spinor",
        max_memory=max_memory,
        ioblk_size=ioblk_size,
        aosym="s1",
    )
    r_outcore.general(
        mol,
        (mo_coeff_L, mo_coeff_S, mo_coeff_S, mo_coeff_L),
        breit_LSSL % prefix,
        intor="int2e_breit_ssp1sps2_spinor",
        max_memory=max_memory,
        ioblk_size=ioblk_size,
        aosym="s1",
    )
    r_outcore.general(
        mol,
        (mo_coeff_S, mo_coeff_L, mo_coeff_L, mo_coeff_S),
        breit_SLLS % prefix,
        intor="int2e_breit_sps1ssp2_spinor",
        max_memory=max_memory,
        ioblk_size=ioblk_size,
        aosym="s1",
    )


def _r_outcore_Gaunt(
    mol, my_RDHF, npes, prefix, max_memory=MAX_MEMORY, ioblk_size=IOBLK_SIZE
):

    from pyscf.ao2mo import r_outcore

    n2c = mol.nao_2c()
    mo_coeff = my_RDHF.mo_coeff
    mo_coeff_mat = numpy.matrix(mo_coeff)

    mo_coeff_pes = mo_coeff_mat[:, n2c : n2c + npes]
    mo_coeff_L = mo_coeff_pes[:n2c, :]
    mo_coeff_S = mo_coeff_pes[n2c:, :]

    r_outcore.general(
        mol,
        (mo_coeff_L, mo_coeff_S, mo_coeff_L, mo_coeff_S),
        breit_LSLS % prefix,
        intor="int2e_ssp1ssp2_spinor",
        max_memory=max_memory,
        ioblk_size=ioblk_size,
        aosym="s1",
    )
    r_outcore.general(
        mol,
        (mo_coeff_S, mo_coeff_L, mo_coeff_S, mo_coeff_L),
        breit_SLSL % prefix,
        intor="int2e_sps1sps2_spinor",
        max_memory=max_memory,
        ioblk_size=ioblk_size,
        aosym="s1",
    )
    r_outcore.general(
        mol,
        (mo_coeff_L, mo_coeff_S, mo_coeff_S, mo_coeff_L),
        breit_LSSL % prefix,
        intor="int2e_ssp1sps2_spinor",
        max_memory=max_memory,
        ioblk_size=ioblk_size,
        aosym="s1",
    )
    r_outcore.general(
        mol,
        (mo_coeff_S, mo_coeff_L, mo_coeff_L, mo_coeff_S),
        breit_SLLS % prefix,
        intor="int2e_sps1ssp2_spinor",
        max_memory=max_memory,
        ioblk_size=ioblk_size,
        aosym="s1",
    )


def _dump_2e_outcore(
    fout, n2c, prefix, with_breit, with_gaunt, IsComplex, symmetry="s1", tol=1e-8
):

    import h5py

    feri_coulomb_LLLL = h5py.File(coulomb_LLLL % prefix, "r")
    feri_coulomb_LLSS = h5py.File(coulomb_LLSS % prefix, "r")
    feri_coulomb_SSLL = h5py.File(coulomb_SSLL % prefix, "r")
    feri_coulomb_SSSS = h5py.File(coulomb_SSSS % prefix, "r")

    if with_breit or with_gaunt:
        feri_breit_LSLS = h5py.File(breit_LSLS % prefix, "r")
        feri_breit_SLSL = h5py.File(breit_SLSL % prefix, "r")
        feri_breit_LSSL = h5py.File(breit_LSSL % prefix, "r")
        feri_breit_SLLS = h5py.File(breit_SLLS % prefix, "r")
    else:
        feri_breit_LSLS = None
        feri_breit_SLSL = None
        feri_breit_LSSL = None
        feri_breit_SLLS = None

    c1 = 0.5 / lib.param.LIGHT_SPEED

    if with_breit:
        _sign_ = 1.0
    else:
        if with_gaunt:
            _sign_ = -1.0

    if symmetry == "s1":
        if IsComplex:
            for i in range(n2c):
                for j in range(n2c):

                    ij = i * n2c + j

                    eri_coulomb = numpy.array(feri_coulomb_LLLL["eri_mo"][ij]).reshape(
                        n2c, n2c
                    )
                    eri_coulomb += (
                        numpy.array(feri_coulomb_LLSS["eri_mo"][ij]).reshape(n2c, n2c)
                        * c1**2
                    )
                    eri_coulomb += (
                        numpy.array(feri_coulomb_SSLL["eri_mo"][ij]).reshape(n2c, n2c)
                        * c1**2
                    )
                    eri_coulomb += (
                        numpy.array(feri_coulomb_SSSS["eri_mo"][ij]).reshape(n2c, n2c)
                        * c1**4
                    )

                    # max_coulomb = numpy.max(numpy.abs(eri_coulomb))

                    # if max_coulomb > tol:

                    for k in range(n2c):
                        # max_coulomb = numpy.max(numpy.abs(eri_coulomb[k]))
                        # if max_coulomb > tol:
                        for l in range(n2c):
                            if abs(eri_coulomb[k][l]) > tol:
                                fout.write(
                                    "%18.12E %18.12E %4d %4d %4d %4d\n"
                                    % (
                                        eri_coulomb[k][l].real,
                                        eri_coulomb[k][l].imag,
                                        i + 1,
                                        j + 1,
                                        k + 1,
                                        l + 1,
                                    )
                                )

                    if with_breit:

                        eri_breit = (
                            numpy.array(feri_breit_LSLS["eri_mo"][ij]).reshape(n2c, n2c)
                            * c1**2
                        )
                        eri_breit += (
                            numpy.array(feri_breit_SLSL["eri_mo"][ij]).reshape(n2c, n2c)
                            * c1**2
                        )
                        eri_breit += (
                            numpy.array(feri_breit_LSSL["eri_mo"][ij]).reshape(n2c, n2c)
                            * c1**2
                        )
                        eri_breit += (
                            numpy.array(feri_breit_SLLS["eri_mo"][ij]).reshape(n2c, n2c)
                            * c1**2
                        )

                        # max_breit = numpy.max(numpy.abs(eri_breit))

                        eri_breit *= _sign_

                        # if max_breit > tol:
                        for k in range(n2c):
                            # max_breit = numpy.max(numpy.abs(eri_breit[k]))
                            # if max_breit > tol:
                            for l in range(n2c):
                                if abs(eri_breit[k][l]) > tol:
                                    fout.write(
                                        "%18.12E %18.12E %4d %4d %4d %4d\n"
                                        % (
                                            eri_breit[k][l].real,
                                            eri_breit[k][l].imag,
                                            n2c + i + 1,
                                            n2c + j + 1,
                                            n2c + k + 1,
                                            n2c + l + 1,
                                        )
                                    )

        else:

            for i in range(n2c):
                for j in range(n2c):

                    ij = i * n2c + j

                    eri_coulomb = numpy.array(feri_coulomb_LLLL["eri_mo"][ij]).reshape(
                        n2c, n2c
                    )
                    eri_coulomb += (
                        numpy.array(feri_coulomb_LLSS["eri_mo"][ij]).reshape(n2c, n2c)
                        * c1**2
                    )
                    eri_coulomb += (
                        numpy.array(feri_coulomb_SSLL["eri_mo"][ij]).reshape(n2c, n2c)
                        * c1**2
                    )
                    eri_coulomb += (
                        numpy.array(feri_coulomb_SSSS["eri_mo"][ij]).reshape(n2c, n2c)
                        * c1**4
                    )

                    # max_coulomb = numpy.max(numpy.abs(eri_coulomb))

                    # if max_coulomb > tol:
                    for k in range(n2c):
                        # max_coulomb = numpy.max(numpy.abs(eri_coulomb[k]))
                        # if max_coulomb > tol:
                        for l in range(n2c):
                            if abs(eri_coulomb[k][l]) > tol:
                                fout.write(
                                    "%18.12E %4d %4d %4d %4d\n"
                                    % (
                                        eri_coulomb[k][l].real,
                                        i + 1,
                                        j + 1,
                                        k + 1,
                                        l + 1,
                                    )
                                )

                    if with_breit:
                        eri_breit = (
                            numpy.array(feri_breit_LSLS["eri_mo"][ij]).reshape(n2c, n2c)
                            * c1**2
                        )
                        eri_breit += (
                            numpy.array(feri_breit_SLSL["eri_mo"][ij]).reshape(n2c, n2c)
                            * c1**2
                        )
                        eri_breit += (
                            numpy.array(feri_breit_LSSL["eri_mo"][ij]).reshape(n2c, n2c)
                            * c1**2
                        )
                        eri_breit += (
                            numpy.array(feri_breit_SLLS["eri_mo"][ij]).reshape(n2c, n2c)
                            * c1**2
                        )

                        # max_breit = numpy.max(numpy.abs(eri_breit))

                        eri_breit *= _sign_

                        # if max_breit > tol:
                        for k in range(n2c):
                            # max_breit = numpy.max(numpy.abs(eri_breit[k]))
                            # if max_breit > tol:
                            for l in range(n2c):
                                if abs(eri_breit[k][l]) > tol:
                                    fout.write(
                                        "%18.12E %4d %4d %4d %4d\n"
                                        % (
                                            eri_breit[k][l].real,
                                            n2c + i + 1,
                                            n2c + j + 1,
                                            n2c + k + 1,
                                            n2c + l + 1,
                                        )
                                    )

    elif symmetry == "s4":

        if IsComplex:
            for i in range(n2c):
                for j in range(i + 1):

                    ij = i * n2c + j

                    eri_coulomb = numpy.array(feri_coulomb_LLLL["eri_mo"][ij]).reshape(
                        n2c, n2c
                    )
                    eri_coulomb += (
                        numpy.array(feri_coulomb_LLSS["eri_mo"][ij]).reshape(n2c, n2c)
                        * c1**2
                    )
                    eri_coulomb += (
                        numpy.array(feri_coulomb_SSLL["eri_mo"][ij]).reshape(n2c, n2c)
                        * c1**2
                    )
                    eri_coulomb += (
                        numpy.array(feri_coulomb_SSSS["eri_mo"][ij]).reshape(n2c, n2c)
                        * c1**4
                    )

                    # max_coulomb = numpy.max(numpy.abs(eri_coulomb))

                    # if max_coulomb > tol:
                    for k in range(i + 1):
                        # max_coulomb = numpy.max(numpy.abs(eri_coulomb[k]))
                        # if max_coulomb > tol:
                        for l in range(n2c):
                            if abs(eri_coulomb[k][l]) > tol:
                                fout.write(
                                    "%18.12E %18.12E %4d %4d %4d %4d\n"
                                    % (
                                        eri_coulomb[k][l].real,
                                        eri_coulomb[k][l].imag,
                                        i + 1,
                                        j + 1,
                                        k + 1,
                                        l + 1,
                                    )
                                )

                    if with_breit:

                        eri_breit = (
                            numpy.array(feri_breit_LSLS["eri_mo"][ij]).reshape(n2c, n2c)
                            * c1**2
                        )
                        eri_breit += (
                            numpy.array(feri_breit_SLSL["eri_mo"][ij]).reshape(n2c, n2c)
                            * c1**2
                        )
                        eri_breit += (
                            numpy.array(feri_breit_LSSL["eri_mo"][ij]).reshape(n2c, n2c)
                            * c1**2
                        )
                        eri_breit += (
                            numpy.array(feri_breit_SLLS["eri_mo"][ij]).reshape(n2c, n2c)
                            * c1**2
                        )

                        # max_breit = numpy.max(numpy.abs(eri_breit))

                        eri_breit *= _sign_

                        # if max_breit > tol:
                        for k in range(i + 1):
                            # max_breit = numpy.max(numpy.abs(eri_breit[k]))
                            # if max_breit > tol:
                            for l in range(n2c):
                                if abs(eri_breit[k][l]) > tol:
                                    fout.write(
                                        "%18.12E %18.12E %4d %4d %4d %4d\n"
                                        % (
                                            eri_breit[k][l].real,
                                            eri_breit[k][l].imag,
                                            n2c + i + 1,
                                            n2c + j + 1,
                                            n2c + k + 1,
                                            n2c + l + 1,
                                        )
                                    )

        else:
            for i in range(n2c):
                for j in range(i + 1):

                    ij = i * n2c + j

                    eri_coulomb = numpy.array(feri_coulomb_LLLL["eri_mo"][ij]).reshape(
                        n2c, n2c
                    )
                    eri_coulomb += (
                        numpy.array(feri_coulomb_LLSS["eri_mo"][ij]).reshape(n2c, n2c)
                        * c1**2
                    )
                    eri_coulomb += (
                        numpy.array(feri_coulomb_SSLL["eri_mo"][ij]).reshape(n2c, n2c)
                        * c1**2
                    )
                    eri_coulomb += (
                        numpy.array(feri_coulomb_SSSS["eri_mo"][ij]).reshape(n2c, n2c)
                        * c1**4
                    )

                    # max_coulomb = numpy.max(numpy.abs(eri_coulomb))

                    # if max_coulomb > tol:
                    for k in range(i + 1):
                        # max_coulomb = numpy.max(numpy.abs(eri_coulomb[k]))
                        # if max_coulomb > tol:
                        for l in range(n2c):
                            if abs(eri_coulomb[k][l]) > tol:
                                fout.write(
                                    "%18.12E %4d %4d %4d %4d\n"
                                    % (
                                        eri_coulomb[k][l].real,
                                        i + 1,
                                        j + 1,
                                        k + 1,
                                        l + 1,
                                    )
                                )

                    if with_breit:
                        eri_breit = (
                            numpy.array(feri_breit_LSLS["eri_mo"][ij]).reshape(n2c, n2c)
                            * c1**2
                        )
                        eri_breit += (
                            numpy.array(feri_breit_SLSL["eri_mo"][ij]).reshape(n2c, n2c)
                            * c1**2
                        )
                        eri_breit += (
                            numpy.array(feri_breit_LSSL["eri_mo"][ij]).reshape(n2c, n2c)
                            * c1**2
                        )
                        eri_breit += (
                            numpy.array(feri_breit_SLLS["eri_mo"][ij]).reshape(n2c, n2c)
                            * c1**2
                        )

                        # max_breit = numpy.max(numpy.abs(eri_breit))

                        eri_breit *= _sign_

                        # if max_breit > tol:
                        for k in range(i + 1):
                            # max_breit = numpy.max(numpy.abs(eri_breit[k]))
                            # if max_breit > tol:
                            for l in range(n2c):
                                if abs(eri_breit[k][l]) > tol:
                                    fout.write(
                                        "%18.12E %4d %4d %4d %4d\n"
                                        % (
                                            eri_breit[k][l].real,
                                            n2c + i + 1,
                                            n2c + j + 1,
                                            n2c + k + 1,
                                            n2c + l + 1,
                                        )
                                    )

    else:
        raise ValueError("Unknown symmetry %s" % symmetry)

    feri_coulomb_LLLL.close()
    feri_coulomb_LLSS.close()
    feri_coulomb_SSLL.close()
    feri_coulomb_SSSS.close()

    if with_breit:
        feri_breit_LSLS.close()
        feri_breit_SLSL.close()
        feri_breit_LSSL.close()
        feri_breit_SLLS.close()


def _is_paldus_123(i, j, k, l):
    return (i != j) and (i != k) and (i != l) and (j != k) and (j != l) and (k != l)


def _is_paldus_456(i, j, k, l):
    return ((i == k) and (j != l) and (i != j)) or ((j == l) and (i != k) and (i != j))


def _is_dg(i, j, k, l):
    is_true_case1 = (i == j) and (k == l)
    is_true_case2 = (i == l) and (j == k)
    return is_true_case1 or is_true_case2
    # return (i == j) and (i == k) and (i == l)


def _dump_2e(
    fout,
    int2e_coulomb,
    int2e_breit,
    with_breit,
    IsComplex,
    symmetry="s1",
    tol=1e-8,
    retain_paldus_123=True,
    retain_paldus_456=True,
    only_dg=False,
):
    """Dump the 2-electron integrals in FCIDUMP format (**incore** mode)

    Args:
        fout: the file object to dump the integrals
        int2e_coulomb: the 2-electron Coulomb integrals
        int2e_breit: the 2-electron Breit integrals
        with_breit: whether to include Breit term
        IsComplex: whether the integrals are complex
        symmetry: the symmetry of the integrals

    Kwargs:

    Returns:

    """

    # tol = 1e-10
    n2c = int2e_coulomb.shape[0]

    if symmetry == "s1":
        if IsComplex:
            for i in range(n2c):
                for j in range(n2c):
                    for k in range(n2c):
                        for l in range(n2c):

                            is_paldus_123 = _is_paldus_123(
                                i // 2, j // 2, k // 2, l // 2
                            )
                            is_paldus_456 = _is_paldus_456(
                                i // 2, j // 2, k // 2, l // 2
                            )
                            is_dg = _is_dg(i // 2, j // 2, k // 2, l // 2)

                            if (not retain_paldus_123) and is_paldus_123:
                                continue
                            if (not retain_paldus_456) and is_paldus_456:
                                continue
                            if (not is_dg) and only_dg:
                                continue

                            if abs(int2e_coulomb[i][j][k][l]) > tol:
                                fout.write(
                                    "%18.12E %18.12E %4d %4d %4d %4d\n"
                                    % (
                                        int2e_coulomb[i][j][k][l].real,
                                        int2e_coulomb[i][j][k][l].imag,
                                        i + 1,
                                        j + 1,
                                        k + 1,
                                        l + 1,
                                    )
                                )
                            if with_breit:
                                if abs(int2e_breit[i][j][k][l]) > tol:
                                    fout.write(
                                        "%18.12E %18.12E %4d %4d %4d %4d\n"
                                        % (
                                            int2e_breit[i][j][k][l].real,
                                            int2e_breit[i][j][k][l].imag,
                                            n2c + i + 1,
                                            n2c + j + 1,
                                            n2c + k + 1,
                                            n2c + l + 1,
                                        )
                                    )
        else:
            for i in range(n2c):
                for j in range(n2c):
                    for k in range(n2c):
                        for l in range(n2c):

                            is_paldus_123 = _is_paldus_123(
                                i // 2, j // 2, k // 2, l // 2
                            )
                            is_paldus_456 = _is_paldus_456(
                                i // 2, j // 2, k // 2, l // 2
                            )
                            is_dg = _is_dg(i // 2, j // 2, k // 2, l // 2)

                            if (not retain_paldus_123) and is_paldus_123:
                                continue
                            if (not retain_paldus_456) and is_paldus_456:
                                continue
                            if (not is_dg) and only_dg:
                                continue

                            if abs(int2e_coulomb[i][j][k][l]) > tol:
                                fout.write(
                                    "%18.12E %4d %4d %4d %4d\n"
                                    % (
                                        int2e_coulomb[i][j][k][l].real,
                                        i + 1,
                                        j + 1,
                                        k + 1,
                                        l + 1,
                                    )
                                )
                            if with_breit:
                                if abs(int2e_breit[i][j][k][l]) > tol:
                                    fout.write(
                                        "%18.12E %4d %4d %4d %4d\n"
                                        % (
                                            int2e_breit[i][j][k][l].real,
                                            n2c + i + 1,
                                            n2c + j + 1,
                                            n2c + k + 1,
                                            n2c + l + 1,
                                        )
                                    )

    elif symmetry == "s4":

        if IsComplex:
            for i in range(n2c):
                for j in range(i + 1):
                    for k in range(i + 1):
                        for l in range(n2c):

                            is_paldus_123 = _is_paldus_123(
                                i // 2, j // 2, k // 2, l // 2
                            )
                            is_paldus_456 = _is_paldus_456(
                                i // 2, j // 2, k // 2, l // 2
                            )
                            is_dg = _is_dg(i // 2, j // 2, k // 2, l // 2)

                            if (not retain_paldus_123) and is_paldus_123:
                                continue
                            if (not retain_paldus_456) and is_paldus_456:
                                continue
                            if (not is_dg) and only_dg:
                                continue

                            if abs(int2e_coulomb[i][j][k][l]) > tol:
                                fout.write(
                                    "%18.12E %18.12E %4d %4d %4d %4d\n"
                                    % (
                                        int2e_coulomb[i][j][k][l].real,
                                        int2e_coulomb[i][j][k][l].imag,
                                        i + 1,
                                        j + 1,
                                        k + 1,
                                        l + 1,
                                    )
                                )
                            if with_breit:
                                if abs(int2e_breit[i][j][k][l]) > tol:
                                    fout.write(
                                        "%18.12E %18.12E %4d %4d %4d %4d\n"
                                        % (
                                            int2e_breit[i][j][k][l].real,
                                            int2e_breit[i][j][k][l].imag,
                                            n2c + i + 1,
                                            n2c + j + 1,
                                            n2c + k + 1,
                                            n2c + l + 1,
                                        )
                                    )
        else:
            for i in range(n2c):
                for j in range(i + 1):
                    for k in range(i + 1):
                        for l in range(n2c):

                            is_paldus_123 = _is_paldus_123(
                                i // 2, j // 2, k // 2, l // 2
                            )
                            is_paldus_456 = _is_paldus_456(
                                i // 2, j // 2, k // 2, l // 2
                            )
                            is_dg = _is_dg(i // 2, j // 2, k // 2, l // 2)

                            if (not retain_paldus_123) and is_paldus_123:
                                continue
                            if (not retain_paldus_456) and is_paldus_456:
                                continue
                            if (not is_dg) and only_dg:
                                continue

                            if abs(int2e_coulomb[i][j][k][l]) > tol:
                                fout.write(
                                    "%18.12E %4d %4d %4d %4d\n"
                                    % (
                                        int2e_coulomb[i][j][k][l].real,
                                        i + 1,
                                        j + 1,
                                        k + 1,
                                        l + 1,
                                    )
                                )
                            if with_breit:
                                if abs(int2e_breit[i][j][k][l]) > tol:
                                    fout.write(
                                        "%18.12E %4d %4d %4d %4d\n"
                                        % (
                                            int2e_breit[i][j][k][l].real,
                                            n2c + i + 1,
                                            n2c + j + 1,
                                            n2c + k + 1,
                                            n2c + l + 1,
                                        )
                                    )

    else:
        raise ValueError("Unknown symmetry %s" % symmetry)
