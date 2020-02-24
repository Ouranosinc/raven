import pytest
import datetime as dt
import numpy as np
import xarray as xr
import netCDF4 as nc

from pywps import Service
from pywps.tests import assert_response_success

from . common import client_for, TESTDATA, CFG_FILE, get_output, urlretrieve
from raven.processes import RavenHMETSProcess
import pdb

class TestRavenERA5Process:

    def test_simple(self):
        client = client_for(Service(processes=[RavenHMETSProcess(), ], cfgfiles=CFG_FILE))

        startYear=2007
        
        params = '9.5019, 0.2774, 6.3942, 0.6884, 1.2875, 5.4134, 2.3641, 0.0973, 0.0464, 0.1998, 0.0222, -1.0919, ' \
                 '2.6851, 0.3740, 1.0000, 0.4739, 0.0114, 0.0243, 0.0069, 310.7211, 916.1947'
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
            
            # Now we need to stack years one after the other
        
        
        main=tmaxYear.merge(tminYear,compat='override')
        main=main.merge(prYear,compat='override')
        main.to_netcdf('/home/ets/TEST_NRCAN.nc')
        pdb.set_trace()
        D = nc.Dataset('/home/ets/TEST_NRCAN.nc', "a")
        D.variables["time"].units = "days since " + str(startYear) + "-01-01 00:00:00"
        D.close()     
        D = nc.Dataset('/home/ets/TEST_NRCAN.nc', "r+")
        D.variables["time"] = list(range(0,D.variables["tasmax"].shape[0]))
        D.close()     
        
        '''
        need to write netcdf file here, I did it outside Python. Probably possible to pass a netcdf file directly?
        '''
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
            .format(ts='/home/ets/TEST_NRCAN.nc', # This is a file on disk, need to pass 'g'
                    params=params,
                    start_date=dt.datetime(2008, 1, 1),
                    end_date=dt.datetime(2008, 8, 10),
                    init='155,455',
                    name='Salmon',
                    run_name='test-hmets-NRCAN',
                    area='4250.6',
                    elevation='843.0',
                    latitude=salmon_lat,
                    longitude=salmon_lon,
                    )
        
        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='raven-hmets',
            datainputs=datainputs)
        pdb.set_trace()
        assert_response_success(resp)
        out = get_output(resp.xml)
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
