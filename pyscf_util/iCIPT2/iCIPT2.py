from pyscf_util.misc.icipt2_inputfile_generator import _Generate_InputFile_iCI
import os
from copy import deepcopy

_FILE_NOT_REMOVE = {
    "start_with": ["FCIDUMP"],
    "end_with": [".py", ".out", ".inp"],
}

from contextlib import contextmanager


@contextmanager
def temporary_update(dictionary, key, value):
    original_value = deepcopy(dictionary.get(key))
    if isinstance(value, str):
        value = [value]
    if value is not None:
        dictionary[key].extend(value)
    try:
        yield
    finally:
        if original_value is None:
            del dictionary[key]
        else:
            dictionary[key] = original_value


class _iCIPT2_Driver:
    def __init__(self, APP, file_not_remove=_FILE_NOT_REMOVE):
        self._APP = APP
        self._file_not_remove = file_not_remove

    def run(
        self,
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
        inputfilename = task_name + ".inp"
        outputfilename = task_name + ".out"
        errorfilename = task_name + ".err"

        if fcidump is None:
            fcidump = "FCIDUMP"

        # (1) generate input file #

        _Generate_InputFile_iCI(
            inputfilename,
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
        )

        # (2) run the application #

        ret_val = os.system(
            "%s %s %s 1>%s 2>%s"
            % (self._APP, inputfilename, fcidump, outputfilename, errorfilename)
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


ICIPT2_CSF_DRIVER = _iCIPT2_Driver(_load_app("ICI_CSF_CPP"), _FILE_NOT_REMOVE)
ICIPT2_DET_DRIVER = _iCIPT2_Driver(_load_app("ICI_DET_CPP"), _FILE_NOT_REMOVE)


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
        ICIPT2_CSF_DRIVER.run(
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
        ICIPT2_DET_DRIVER.run(
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
