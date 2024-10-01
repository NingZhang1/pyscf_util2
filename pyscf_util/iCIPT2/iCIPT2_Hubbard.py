from pyscf_util.iCIPT2.iCIPT2 import temporary_update, _FILE_NOT_REMOVE
import os
from pyscf_util.misc.icipt2_inputfile_generator import _Generate_InputFile_iCI


class _iCIPT2_Hubbard_Driver:
    def __init__(self, APP, dimension, file_not_remove=_FILE_NOT_REMOVE):
        self._APP = APP
        self._dimension = dimension
        self._file_not_remove = file_not_remove

    def run(
        self,
        U,
        Lx,
        Ly,
        task_name,
        segment,
        nelec_val,
        rotatemo=0,
        cmin: str | float = 1e-4,
        perturbation=0,
        Task: str = None,
        etol=1e-7,
        selection=1,
        direct=None,
        start_with=None,
        end_with=None,
    ):
        inputfilename = task_name + ".inp"
        outputfilename = task_name + ".out"
        errorfilename = task_name + ".err"

        # (1) generate input file #

        _Generate_InputFile_iCI(
            inputfilename,
            segment,
            nelec_val,
            rotatemo,
            cmin,
            perturbation,
            0,
            0,
            Task,
            0,
            etol,
            selection,
            0,
            direct,
        )

        # (2) run the application #

        if self._dimension == 1:
            ret_val = os.system(
                "%s %s %f 1>%s 2>%s"
                % (self._APP, inputfilename, U, outputfilename, errorfilename)
            )
        else:
            ret_val = os.system(
                "%s %s %f %d %d 1>%s 2>%s"
                % (self._APP, inputfilename, Lx, Ly, U, outputfilename, errorfilename)
            )

        if ret_val != 0:
            raise ValueError("Error in running the application")

        # (3) remove the input file #

        for file in os.listdir():
            if not os.path.isfile(file):
                continue
            with temporary_update(self._file_not_remove, "start_with", start_with):
                if any(file.endswith(ext) for ext in self._file_not_remove["end_with"]):
                    continue
            with temporary_update(self._file_not_remove, "end_with", end_with):
                if any(
                    file.startswith(start)
                    for start in self._file_not_remove["start_with"]
                ):
                    continue
            os.remove(file)


def _load_app(env_var):
    app = os.getenv(env_var)
    if app is None:
        raise ValueError(f"Environment variable {env_var} is not set")
    return app


# DRIVER #

ICIPT2_Hubbard_CSF_1D_DRIVER = _iCIPT2_Hubbard_Driver(
    _load_app("ICIPT2_Hubbard_CSF_1D_CPP"), 1, _FILE_NOT_REMOVE
)

ICIPT2_Hubbard_CSF_2D_DRIVER = _iCIPT2_Hubbard_Driver(
    _load_app("ICIPT2_Hubbard_CSF_2D_CPP"), 2, _FILE_NOT_REMOVE
)


def kernel_1D(
    U,
    task_name,
    segment,
    nelec_val,
    rotatemo=0,
    cmin: str | float = 1e-4,
    perturbation=0,
    Task: str = None,
    etol=1e-7,
    selection=1,
    direct=None,
    start_with=None,
    end_with=None,
):
    ICIPT2_Hubbard_CSF_1D_DRIVER.run(
        U,
        0,
        0,
        task_name,
        segment,
        nelec_val,
        rotatemo,
        cmin,
        perturbation,
        Task,
        etol,
        selection,
        direct,
        start_with,
        end_with,
    )


def kernel_2D(
    U,
    Lx,
    Ly,
    task_name,
    segment,
    nelec_val,
    rotatemo=0,
    cmin: str | float = 1e-4,
    perturbation=0,
    Task: str = None,
    etol=1e-7,
    selection=1,
    direct=None,
    start_with=None,
    end_with=None,
):
    ICIPT2_Hubbard_CSF_2D_DRIVER.run(
        U,
        Lx,
        Ly,
        task_name,
        segment,
        nelec_val,
        rotatemo,
        cmin,
        perturbation,
        Task,
        etol,
        selection,
        direct,
        start_with,
        end_with,
    )


if __name__ == "__main__":
    kernel_1D(8.0, "test1d", "0 0 0 8 0 0", 8, Task="0 0 1 1", perturbation=1)
    kernel_2D(8.0, 2, 4, "test2d", "0 0 0 8 0 0", 8, Task="0 0 1 1", perturbation=1)
