import pytest
import datetime as dt
import numpy as np
import xarray as xr
import json
import netCDF4 as nc
from pywps import Service
from pywps.tests import assert_response_success
import os

from . common import client_for, TESTDATA, CFG_FILE, get_output, urlretrieve
from raven.processes import RavenHMETSProcess
from raven.models import HMETS



# THIS SECTION BASICALLY PREPARES A NETCDF FILE TO RUN THE TESTS...
startYear=2007
path = os.getcwd()
filepath = path + "/NRCAN_ts.nc"
tasmax='https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/birdhouse/nrcan/nrcan_canada_daily_v2/tasmax/nrcan_canada_daily_tasmax_'
tasmin = 'https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/birdhouse/nrcan/nrcan_canada_daily_v2/tasmin/nrcan_canada_daily_tasmin_'
precip = 'https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/birdhouse/nrcan/nrcan_canada_daily_v2/pr/nrcan_canada_daily_pr_'

# Get information for given catchment, could be passed in parameter to the function
salmon=xr.open_dataset(TESTDATA['raven-hmets-nc-ts'])
salmon_lat = salmon.lat.values[0]
salmon_lon = salmon.lon.values[0]

# Get first year
firstYear=str(startYear)
tmaxYear=xr.open_dataset(tasmax + firstYear + '.nc').sel(lon=salmon_lon, lat=salmon_lat, method='nearest')
tminYear=xr.open_dataset(tasmin + firstYear + '.nc').sel(lon=salmon_lon, lat=salmon_lat, method='nearest')
prYear=xr.open_dataset(precip + firstYear + '.nc').sel(lon=salmon_lon, lat=salmon_lat, method='nearest')

# Now extract all years
for i in range(startYear+1,2009):

    tmaxYear=xr.concat([tmaxYear,xr.open_dataset(tasmax + str(i) + '.nc').sel(lon=salmon_lon, lat=salmon_lat, method='nearest')],dim='time')
    tminYear=xr.concat([tminYear,xr.open_dataset(tasmin + str(i) + '.nc').sel(lon=salmon_lon, lat=salmon_lat, method='nearest')],dim='time')
    prYear=xr.concat([prYear,xr.open_dataset(precip + str(i) + '.nc').sel(lon=salmon_lon, lat=salmon_lat, method='nearest')],dim='time')

# Merge and fix time units so that RAVEN understands
main=tmaxYear.merge(tminYear,compat='override')
main=main.merge(prYear,compat='override')
main.to_netcdf(filepath)
D = nc.Dataset(filepath, "a")
D.variables["time"].units = "days since " + str(startYear) + "-01-01 00:00:00"
D.close()     
D = nc.Dataset(filepath, "r+")
D.variables["time"] = list(range(0,D.variables["tasmax"].shape[0]))
D.close() 


class TestRavenNRCAN:
    def test_simple(self):

        
        params = (9.5019, 0.2774, 6.3942, 0.6884, 1.2875, 5.4134, 2.3641, 0.0973, 0.0464, 0.1998, 0.0222, -1.0919,
          2.6851, 0.3740, 1.0000, 0.4739, 0.0114, 0.0243, 0.0069, 310.7211, 916.1947)
                 
        model = HMETS()
        model(ts=filepath,
              params=params,
              start_date=dt.datetime(2007, 1, 1),
              end_date=dt.datetime(2007, 8, 10),
              name='Salmon',
              run_name='test-hmets-NRCAN',
              area=4250.6,
              elevation=843.0,
              latitude=54.4848,
              longitude=-123.3659,
              rain_snow_fraction="RAINSNOW_DINGMAN",
              tasmax={'linear_transform': (1.0, -273.15)},
              tasmin={'linear_transform': (1.0, -273.15)},
              pr={'linear_transform': (86400, 0.0)}
              )
        
       
        
class TestRavenERA5Process:

    def test_simple(self):
        client = client_for(Service(processes=[RavenHMETSProcess(), ], cfgfiles=CFG_FILE))

        params = '9.5019, 0.2774, 6.3942, 0.6884, 1.2875, 5.4134, 2.3641, 0.0973, 0.0464, 0.1998, 0.0222, -1.0919, ' \
                 '2.6851, 0.3740, 1.0000, 0.4739, 0.0114, 0.0243, 0.0069, 310.7211, 916.1947'
        
        datainputs = "ts=files@xlink:href=file://{ts};" \
                     "params={params};" \
                     "start_date={start_date};" \
                     "end_date={end_date};" \
                     "init={init};" \
                     "name={name};" \
                     "run_name={run_name};" \
                     "area={area};" \
                     "latitude={latitude};" \
                     "longitude={longitude};" \
                     "elevation={elevation};" \
                     "rain_snow_fraction={rain_snow_fraction};"\
                     "nc_spec={tasmax};"\
                     "nc_spec={tasmin};"\
                     "nc_spec={pr}"\
            .format(ts=filepath, # This is a file on disk, need to pass 'g'
                    params=params,
                    start_date=dt.datetime(2007, 1, 1),
                    end_date=dt.datetime(2007, 8, 10),
                    init='155,455',
                    name='Salmon',
                    run_name='test-hmets-NRCAN',
                    area='4250.6',
                    elevation='843.0',
                    latitude=salmon_lat,
                    longitude=salmon_lon,
                    rain_snow_fraction="RAINSNOW_DINGMAN",
                    tasmax=json.dumps({'tasmax': {'linear_transform': (1.0, -273.15)}}),
                    tasmin=json.dumps({'tasmin': {'linear_transform': (1.0, -273.15)}}),
                    pr=json.dumps({'pr': {'linear_transform': (86400.0, 0.0)}}),
                    )
            
        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='raven-hmets',
            datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)
        
        """
        # THERE IS NO QOBS HERE SO NO TESTING OF NASH PERFORMANCE
        
        assert 'diagnostics' in out
        tmp_file, _ = urlretrieve(out['diagnostics'])
        tmp_content = open(tmp_file).readlines()
        

        # checking correctness of NSE (full period 1954-2011 would be NSE=0.636015 as template in Wiki)
        assert 'DIAG_NASH_SUTCLIFFE' in tmp_content[0]
        idx_diag = tmp_content[0].split(',').index("DIAG_NASH_SUTCLIFFE")
        diag = np.float(tmp_content[1].split(',')[idx_diag])
        np.testing.assert_almost_equal(diag, -7.03141, 4, err_msg='NSE is not matching expected value')

        # checking correctness of RMSE (full period 1954-2011 would be RMSE=28.3759 as template in wiki)
        assert 'DIAG_RMSE' in tmp_content[0]
        idx_diag = tmp_content[0].split(',').index("DIAG_RMSE")
        diag = np.float(tmp_content[1].split(',')[idx_diag])

        np.testing.assert_almost_equal(diag, 101.745, 4, err_msg='RMSE is not matching expected value')
        """