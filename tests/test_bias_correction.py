import datetime as dt
import xarray as xr
import xclim.sdba as sdba


class TestBiasCorrect:
    def test_bias_correction(self):
    
        ref_data='https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/datasets/simulations/cmip5/atmos/day_MPI-ESM-LR_historical.ncml'
        fut_data='https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/datasets/simulations/cmip5/atmos/day_MPI-ESM-LR_historical+rcp85.ncml'
        hist_data="https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/birdhouse/1-Datasets/gridded_obs/nrcan_v2.ncml"
                  
        lat=54.484
        lon=-123.36   
        
        # CAREFUL! ERA5 IS NOT ThE SAME LONGITUDE
        # Also, latitude goes from high to low, so I need to have it go from lat+1 to lat-1 in the slice.
        # For the NRCan dataset, I cannot have more than about 10 years as I get a "NetCDF: DAP failure" which I think is related to a timeout.
        ds = (xr.open_dataset(hist_data).sel(lat=slice(lat + 1, lat - 1),lon=slice(lon - 1, lon + 1), time=slice(dt.datetime(1991,1,1), dt.datetime(2010,12,31))).mean(dim={"lat", "lon"}, keep_attrs=True))

        # For lon in 0-360 format, need to add an auto-checker.
        lon = 260
        ds2 = (xr.open_dataset(ref_data).sel(lat=slice(lat - 1, lat + 1),lon=slice(lon - 1, lon + 1), time=slice(dt.datetime(1981,1,1), dt.datetime(2010,12,31))).mean(dim={"lat", "lon"}, keep_attrs=True))
        ds3 = (xr.open_dataset(fut_data).sel(lat=slice(lat - 1, lat + 1),lon=slice(lon - 1, lon + 1), time=slice(dt.datetime(2041,1,1), dt.datetime(2070,12,31))).mean(dim={"lat", "lon"}, keep_attrs=True))
     
        # Here data in ds, ds2 and ds3 are NaN!
        
        group_month_nowindow = sdba.utils.Grouper('time.month')
        Adj = sdba.DetrendedQuantileMapping(nquantiles=50, kind='+', group=group_month_nowindow)
        Adj.train(ds['pr'],ds2['pr'])
        Scen = Adj.adjust(ds3['pr'], interp="linear")
        Adj.ds.af  # adjustment factors.
        
        print(Adj.ds.af)
