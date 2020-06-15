import datetime as dt
import xarray as xr
import xclim.sdba as sdba
from .common import TESTDATA
import pdb

def test_bias_correction():
    
    ref_data='https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/datasets/simulations/cmip5/atmos/day_MPI-ESM-LR_historical.ncml'
    fut_data='https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/datasets/simulations/cmip5/atmos/day_MPI-ESM-LR_historical+rcp85.ncml'
    hist_data='https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/datasets/reanalyses/era5.ncml'
    
    lat=45
    lon=-70
    
    # % CAREFUL! ERA5 IS NOT ThE SAME LONGITUDE!
    ds = (xr.open_dataset(hist_data).sel(lat=slice(lat + 1, lat - 1),lon=slice(lon - 1, lon + 1), time=slice(dt.datetime(1981,1,1), dt.datetime(2010,12,31))).mean(dim={"lat", "lon"}, keep_attrs=True))
    ds.to_netcdf('hist.nc')
    ds2 = (xr.open_dataset(ref_data).sel(lat=slice(lat + 1, lat - 1),lon=slice(lon - 1, lon + 1), time=slice(dt.datetime(1981,1,1), dt.datetime(2010,12,31))).mean(dim={"lat", "lon"}, keep_attrs=True))
    ds2.to_netcdf('ref.nc')
    ds3 = (xr.open_dataset(fut_data).sel(lat=slice(lat + 1, lat - 1),lon=slice(lon - 1, lon + 1), time=slice(dt.datetime(2041,1,1), dt.datetime(2070,12,31))).mean(dim={"lat", "lon"}, keep_attrs=True))
    ds3.to_netcdf('fut.nc')
    
    pdb.set_trace()
    Adj = sdba.adjustment(group="time.month")
    Adj.train('ref.nc','hist.nc')
    scen = Adj.adjust(fut, interp="linear")
    Adj.ds.af  # adjustment factors.
    
   

#    assert (qsim.max() > 1)
#    assert (len(ens) == 2)
#    assert 'realization' in ens.dims
#    assert 'param' in ens.dims
