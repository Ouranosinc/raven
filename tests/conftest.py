import pytest
import xarray as xr
from xclim.indicators.land._streamflow import stats, fit
from raven.tutorial import get_file
from raven.models import Raven
from .common import TESTDATA

SALMON_coords = (-123.3659, 54.4848)  # (lon, lat)
RAVEN = Raven


@pytest.fixture
def q_sim_1(tmp_path):
    """A file storing a Raven streamflow simulation over one basin."""

    return TESTDATA['simfile_single']


@pytest.fixture
def ts_stats(q_sim_1, tmp_path):
    q = xr.open_dataset(q_sim_1).q_sim
    ts = stats(q, op='max')
    fn = tmp_path / 'ts_stats.nc'
    ts.to_netcdf(fn)
    return fn


@pytest.fixture
def params(ts_stats, tmp_path):
    ds = xr.open_dataset(ts_stats)
    name = list(ds.data_vars.keys()).pop()
    q = ds[name]
    p = fit(q, dist='gumbel_r')
    fn = tmp_path / 'fit.nc'
    p.to_netcdf(fn)
    return fn


@pytest.fixture
def era5_hr():
    """Return a netCDF file with hourly ERA5 data at the Salmon location.

    url = "http://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/datasets/reanalyses/era5.ncml"
    ds = xr.open_dataset(url)
    lon, lat = SALMON_coords
    out = ds.sel(longitude=lon + 360, latitude=lat, method='nearest').sel(time=slice("2018-01-01", "2018-1-8"))
    out.to_netcdf("raven-testdata/era5/tas_pr_20180101-20180108.nc")
    """
    return get_file("era5/tas_pr_20180101-20180108.nc")
