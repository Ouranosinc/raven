import pytest
import datetime as dt
import numpy as np
import xarray as xr

from pywps import Service
from pywps.tests import assert_response_success

from . common import client_for, TESTDATA, CFG_FILE, get_output, urlretrieve
from raven.processes import RavenHMETSProcess
import pdb

class TestRavenERA5Process:

    def test_simple(self):
        client = client_for(Service(processes=[RavenHMETSProcess(), ], cfgfiles=CFG_FILE))

        params = '9.5019, 0.2774, 6.3942, 0.6884, 1.2875, 5.4134, 2.3641, 0.0973, 0.0464, 0.1998, 0.0222, -1.0919, ' \
                 '2.6851, 0.3740, 1.0000, 0.4739, 0.0114, 0.0243, 0.0069, 310.7211, 916.1947'
        
        tas = "https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/birdhouse/ecmwf/era5/tas_era5_reanalysis_hourly_2018.nc"
        precip = "https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/birdhouse/ecmwf/era5/pr_era5_reanalysis_hourly_2018.nc"
        tdata=xr.open_dataset(tas)
        pdata=xr.open_dataset(precip)

        salmon=xr.open_dataset(TESTDATA['raven-hmets-nc-ts'])
        salmon_lat = salmon.lat.values[0]
        salmon_lon = salmon.lon.values[0]
        pdb.set_trace()
        
        taslatlon=tdata.sel(longitude=salmon_lon, latitude=salmon_lat, method='nearest')
        prlatlon=pdata.sel(longitude=salmon_lon, latitude=salmon_lat, method='nearest')
        g=taslatlon.merge(prlatlon)
             
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
            .format(ts="/home/ets/testERA5.nc", # This is a file on disk, need to pass 'g'
                    params=params,
                    start_date=dt.datetime(2018, 1, 1),
                    end_date=dt.datetime(2018, 8, 10),
                    init='155,455',
                    name='Salmon',
                    run_name='test-hmets-era5',
                    area='4250.6',
                    elevation='843.0',
                    latitude=54.4848,
                    longitude=-123.3659,
                    )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='raven-hmets',
            datainputs=datainputs)

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
