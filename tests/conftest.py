import pytest
import xarray as xr
from xclim.indicators.land._streamflow import fit, stats

from .common import TESTDATA


@pytest.fixture
def q_sim_1(tmp_path):
    """A file storing a Raven streamflow simulation over one basin."""

    return TESTDATA["simfile_single"]


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
