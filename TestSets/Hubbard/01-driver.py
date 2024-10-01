from pyscf_util.iCIPT2.iCIPT2 import kernel
from pyscf_util.iCIPT2.iCIPT2_Hubbard import kernel_1D, kernel_2D

## first run iCIPT2 with input integrals for Hubbard ##


kernel(
    True,
    task_name="iCIPT2_Hubbard_CSF_1D_16",
    fcidump="FCIDUMP",
    segment="0 2 6 6 2 0",
    nelec_val=12,
    cmin="1e-4 5e-5 3e-5 1.5e-5 9e-6 7e-6 5e-6",
    perturbation=1,
    Task="0 0 1 1",
    inputocfg=0,
)

for L in [16, 20, 24]:
    kernel_1D(
        8.0,
        "test1d_%d" % L,
        "0 0 0 %d 0 0" % (L),
        L,
        Task="0 0 1 1",
        perturbation=1,
        cmin="1e-4 5e-5 3e-5 1.5e-5 9e-6 7e-6 5e-6",
    )

for Lx in [4]:
    for Ly in [4, 5, 6]:
        kernel_2D(
            8.0,
            Lx,
            Ly,
            "test2d_%d_%d" % (Lx, Ly),
            "0 0 0 %d 0 0" % (Lx * Ly),
            (Lx * Ly),
            Task="0 0 1 1",
            perturbation=1,
            cmin="1e-4 5e-5 3e-5 1.5e-5 9e-6 7e-6 5e-6",
        )
