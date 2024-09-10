from pyscf_util.iCIPT2.iCIPT2 import _iCIPT2_Driver, _load_app, _FILE_NOT_REMOVE

ICIPT2_COOV_CSF_DRIVER = _iCIPT2_Driver(_load_app("ICI_CSF_COOV"), _FILE_NOT_REMOVE)
ICIPT2_COOV_DET_DRIVER = _iCIPT2_Driver(_load_app("ICI_DET_COOV"), _FILE_NOT_REMOVE)


def kernel(
    IsCSF: bool,
    task_name,
    fcidump,
    segment,
    nelec_val,
    rotatemo=0,
    cmin: str | float = 1e-4,
    perturbation=0,
    dumprdm=0,
    relative=0,
    Task: str = None,
    inputocfg=0,
    etol=1e-7,
    selection=1,
    doublegroup=None,
    direct=None,
    start_with=None,
    end_with=None,
):
    if IsCSF:
        ICIPT2_COOV_CSF_DRIVER.run(
            task_name,
            fcidump,
            segment,
            nelec_val,
            rotatemo,
            cmin,
            perturbation,
            dumprdm,
            relative,
            Task,
            inputocfg,
            etol,
            selection,
            doublegroup,
            direct,
            start_with,
            end_with,
        )
    else:
        ICIPT2_COOV_DET_DRIVER.run(
            task_name,
            fcidump,
            segment,
            nelec_val,
            rotatemo,
            cmin,
            perturbation,
            dumprdm,
            relative,
            Task,
            inputocfg,
            etol,
            selection,
            doublegroup,
            direct,
            start_with,
            end_with,
        )
