import os
import shutil
from pathlib import Path
from typing import Optional, Union

import pytest
import xarray as xr
from filelock import FileLock
from ravenpy.utilities.testdata import _default_cache_dir
from ravenpy.utilities.testdata import get_file as _get_file
from ravenpy.utilities.testdata import get_local_testdata as _get_local_testdata
from xclim.indicators.land import fit, stats

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

    models = [
        "gr4j-cemaneige",
        "hbvec",
        "hmets",
        "mohyse",
    ]

    data_entries = list()
    entries = [
        "ostrich-{model}/OstRandomNumbers.txt",
        "ostrich-{model}/ostIn.txt",
        "ostrich-{model}/raven-{model0}-salmon.rvc",
        "ostrich-{model}/raven-{model0}-salmon.rvc.tpl",
        "ostrich-{model}/raven-{model0}-salmon.rvh",
        "ostrich-{model}/raven-{model0}-salmon.rvh.tpl",
        "ostrich-{model}/raven-{model0}-salmon.rvi",
        "ostrich-{model}/raven-{model0}-salmon.rvi.tpl",
        "ostrich-{model}/raven-{model0}-salmon.rvp",
        "ostrich-{model}/raven-{model0}-salmon.rvp.tpl",
        "ostrich-{model}/raven-{model0}-salmon.rvt",
        "ostrich-{model}/raven-{model0}-salmon.rvt.tpl",
        "raven-{model}/Salmon-River-Near-Prince-George_Qobs_daily.rvt",
        "raven-{model}/Salmon-River-Near-Prince-George_meteo_daily.rvt",
        "raven-{model}/raven-{model0}-salmon.rvc",
        "raven-{model}/raven-{model0}-salmon.rvh",
        "raven-{model}/raven-{model0}-salmon.rvi",
        "raven-{model}/raven-{model0}-salmon.rvp",
        "raven-{model}/raven-{model0}-salmon.rvt",
    ]
    for model in models:
        for entry in entries:
            data_entries.append(entry.format(model=model, model0=model.split("-")[0]))

    data_entries.extend(
        [
            "XSS_forecast_data/XSS_fcst_det.nc",
            "XSS_forecast_data/XSS_fcst_ens.nc",
            "XSS_forecast_data/XSS_obs.nc",
            "cec_nalcms2010_30m/cec_nalcms_subQC.tiff",
            "donneesqc_mrc_poly/mrc_subset.gml",
            "donneesqc_mrc_poly/mrc_subset.zip",
            "earthenv_dem_90m/earthenv_dem90_southernQuebec.tiff",
            "era5/tas_pr_20180101-20180108.nc",
            "flood_risk/XSS_fcst_det.nc",
            "flood_risk/XSS_fcst_ens.nc",
            "gr4j_cemaneige/solution.rvc",
            "hydro_simulations/raven-gr4j-cemaneige-sim_hmets-0_Hydrographs.nc",
            "nrcan/NRCAN_2006-2007_subset.nc",
            "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc",
            "ts_stats_outputs/out.nc",
            "watershed_vector/Basin_test.zip",
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


@pytest.fixture
def q_sim_1(tmp_path):
    """A file storing a Raven streamflow simulation over one basin."""

    return _get_local_testdata(
        "hydro_simulations/raven-gr4j-cemaneige-sim_hmets-0_Hydrographs.nc",
        temp_folder=tmp_path,
        branch=MAIN_TESTDATA_BRANCH,
        _local_cache=_default_cache_dir,
    )


@pytest.fixture
def ts_stats(q_sim_1, tmp_path):
    q = xr.open_dataset(q_sim_1).q_sim
    ts = stats(q, op="max")
    fn = tmp_path / "ts_stats.nc"
    ts.to_netcdf(fn)
    return fn


@pytest.fixture
def params(ts_stats, tmp_path):
    ds = xr.open_dataset(ts_stats)
    name = list(ds.data_vars.keys()).pop()
    q = ds[name]
    p = fit(q, dist="gumbel_r")
    fn = tmp_path / "fit.nc"
    p.to_netcdf(fn)
    return fn


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
