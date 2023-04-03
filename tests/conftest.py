import os
import shutil
from pathlib import Path
from typing import Optional, Union

import pytest
from filelock import FileLock
from ravenpy.utilities.testdata import _default_cache_dir
from ravenpy.utilities.testdata import get_file as _get_file
from ravenpy.utilities.testdata import get_local_testdata as _get_local_testdata

MAIN_TESTDATA_BRANCH = os.getenv("MAIN_TESTDATA_BRANCH", "master")
SKIP_TEST_DATA = os.getenv("SKIP_TEST_DATA")


def populate_testing_data(
    temp_folder: Optional[Path] = None,
    branch: str = MAIN_TESTDATA_BRANCH,
    _local_cache: Path = _default_cache_dir,
):
    if _local_cache.joinpath(".data_written").exists():
        # This flag prevents multiple calls from re-attempting to download testing data in the same pytest run
        return

    data_entries = list()

    data_entries.extend(
        [
            "cec_nalcms2010_30m/cec_nalcms_subQC.tiff",
            "donneesqc_mrc_poly/mrc_subset.gml",
            "donneesqc_mrc_poly/mrc_subset.zip",
            "earthenv_dem_90m/earthenv_dem90_southernQuebec.tiff",
            "polygons/Basin_10.zip",
            "watershed_vector/Basin_test.zip",
            "watershed_vector/LSJ_LL.zip",
        ]
    )

    data = dict()
    for filepattern in data_entries:
        if temp_folder is None:
            try:
                data[filepattern] = _get_file(
                    filepattern, branch=branch, cache_dir=_local_cache
                )
            except FileNotFoundError:
                continue
        elif temp_folder:
            try:
                data[filepattern] = _get_local_testdata(
                    filepattern,
                    temp_folder=temp_folder,
                    branch=branch,
                    _local_cache=_local_cache,
                )
            except FileNotFoundError:
                continue

    return


@pytest.fixture(scope="session", autouse=True)
def gather_session_data(threadsafe_data_dir, worker_id):
    """Gather testing data on pytest run.
    When running pytest with multiple workers, one worker will copy data remotely to _default_cache_dir while
    other workers wait using lockfile. Once the lock is released, all workers will copy data to their local
    threadsafe_data_dir."""
    if worker_id == "master":
        if not SKIP_TEST_DATA:
            populate_testing_data(branch=MAIN_TESTDATA_BRANCH)
    else:
        if not SKIP_TEST_DATA:
            _default_cache_dir.mkdir(exist_ok=True)
            test_data_being_written = FileLock(_default_cache_dir.joinpath(".lock"))
            with test_data_being_written as fl:
                # This flag prevents multiple calls from re-attempting to download testing data in the same pytest run
                populate_testing_data(branch=MAIN_TESTDATA_BRANCH)
                _default_cache_dir.joinpath(".data_written").touch()
            fl.acquire()
        shutil.copytree(_default_cache_dir, threadsafe_data_dir)


@pytest.fixture(scope="session", autouse=True)
def cleanup(request):
    """Cleanup a testing file once we are finished.
    This flag prevents remote data from being downloaded multiple times in the same pytest run.
    """

    def remove_data_written_flag():
        flag = _default_cache_dir.joinpath(".data_written")
        if flag.exists():
            flag.unlink()

    request.addfinalizer(remove_data_written_flag)


@pytest.fixture(scope="session")
def threadsafe_data_dir(tmp_path_factory) -> Path:
    """Constructor for worker-session temporary data folders."""
    return Path(tmp_path_factory.getbasetemp().joinpath("data"))


@pytest.fixture(scope="session")
def get_file(threadsafe_data_dir):
    def _get_session_scoped_file(file: Union[str, Path]):
        return _get_file(
            file, cache_dir=threadsafe_data_dir, branch=MAIN_TESTDATA_BRANCH
        )

    return _get_session_scoped_file


@pytest.fixture(scope="session")
def get_local_testdata(threadsafe_data_dir):
    def _get_session_scoped_local_testdata(file: Union[str, Path]):
        return _get_local_testdata(
            file,
            temp_folder=threadsafe_data_dir,
            branch=MAIN_TESTDATA_BRANCH,
            _local_cache=_default_cache_dir,
        )

    return _get_session_scoped_local_testdata
