import pytest
import xarray as xr
from urlpath import URL
from xclim.streamflow import stats, fit

from .common import TESTDATA, TD

SALMON_coords = (-123.3659, 54.4848)  # (lon, lat)


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
    import netCDF4 as nc
    """Return a netCDF file with hourly ERA5 data at the Salmon location."""
    path = TD / "ERA5" / "ts.nc"

    if not path.exists():
        # Fetch the data and save to disk if the file has not been created yet.
        path.parent.mkdir(exist_ok=True)
        thredds = URL("https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/birdhouse/ecmwf/era5")
        tas = str(thredds / "tas_era5_reanalysis_hourly_2018.nc")
        pr = str(thredds / "pr_era5_reanalysis_hourly_2018.nc")

        ds = xr.open_mfdataset([tas, pr], combine="by_coords")
        lon, lat = SALMON_coords
        out = ds.sel(longitude=lon + 360, latitude=lat, method='nearest')
        out.to_netcdf(path)

    D = nc.Dataset(path, "a")
    D.variables["time"].units = "hours since 1900-01-01 00:00:00"
    D.close()
    return path
