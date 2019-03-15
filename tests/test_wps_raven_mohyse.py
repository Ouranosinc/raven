import pytest
import datetime as dt
import numpy as np

from pywps import Service
from pywps.tests import assert_response_success

from .common import client_for, TESTDATA, CFG_FILE, get_output, urlretrieve
from raven.processes import RavenMOHYSEProcess
import xarray as xr

class TestRavenMOHYSEProcess:

    def test_simple(self):
        client = client_for(Service(processes=[RavenMOHYSEProcess(), ], cfgfiles=CFG_FILE))

        params = '1.0000, 0.0468, 4.2952, 2.6580, 0.4038, 0.0621, 0.0273, 0.0453'
        hrus = '0.9039, 5.6179775'

        datainputs = "ts=files@xlink:href=file://{ts};" \
                     "params={params};" \
                     "hrus={hrus};" \
                     "start_date={start_date};" \
                     "end_date={end_date};" \
                     "name={name};" \
                     "run_name={run_name};" \
                     "area={area};" \
                     "latitude={latitude};" \
                     "longitude={longitude};" \
                     "elevation={elevation};" \
            .format(ts=TESTDATA['raven-mohyse-nc-ts'],
                    params=params,
                    hrus=hrus,
                    start_date=dt.datetime(2000, 1, 1),
                    end_date=dt.datetime(2002, 1, 1),
                    name='Salmon',
                    run_name='test-mohyse',
                    area='4250.6',
                    elevation='843.0',
                    latitude=54.4848,
                    longitude=-123.3659,
                    )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='raven-mohyse',
            datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)

        assert 'diagnostics' in out
        tmp_file, _ = urlretrieve(out['diagnostics'])
        tmp_content = open(tmp_file).readlines()

        # checking correctness of NSE (full period 1954-2011 would be NSE=0.391103 as template in Wiki)
        assert 'DIAG_NASH_SUTCLIFFE' in tmp_content[0]
        idx_diag = tmp_content[0].split(',').index("DIAG_NASH_SUTCLIFFE")
        diag = np.float(tmp_content[1].split(',')[idx_diag])
        np.testing.assert_almost_equal(diag, 0.194612, 4, err_msg='NSE is not matching expected value')

        # checking correctness of RMSE (full period 1954-2011 would be RMSE=36.7012 as template in wiki)
        assert 'DIAG_RMSE' in tmp_content[0]
        idx_diag = tmp_content[0].split(',').index("DIAG_RMSE")
        diag = np.float(tmp_content[1].split(',')[idx_diag])
        np.testing.assert_almost_equal(diag, 32.2197, 4, err_msg='RMSE is not matching expected value')
        
    def test_parallel(self):
        client = client_for(Service(processes=[RavenMOHYSEProcess(), ], cfgfiles=CFG_FILE))

        params1 = '1.03, 0.046, 4.2952, 2.6580, 0.4038, 0.0621, 0.0273, 0.0453'
        hrus1 = '0.9039, 5.6179775'
        params2 = '1.05, 0.0468, 4.29, 2.65, 0.45, 0.062, 0.027, 0.043'
        hrus2 = '0.903, 5.65'
        params3 = '0.98, 0.04, 4.2, 2.6, 0.40, 0.05, 0.03, 0.03'
        hrus3 = '0.90, 5.5'
        
        datainputs = "ts=files@xlink:href=file://{ts};" \
                     "params={params1};" \
                     "params={params2};" \
                     "params={params3};" \
                     "hrus={hrus1};" \
                     "hrus={hrus2};" \
                     "hrus={hrus3};" \
                     "start_date={start_date};" \
                     "end_date={end_date};" \
                     "name={name};" \
                     "area={area};" \
                     "latitude={latitude};" \
                     "longitude={longitude};" \
                     "elevation={elevation};" \
            .format(ts=TESTDATA['raven-mohyse-nc-ts'],
                    params1=params1,
                    params2=params2,
                    params3=params3,
                    hrus1=hrus1,
                    hrus2=hrus2,
                    hrus3=hrus3,
                    start_date=dt.datetime(2000, 1, 1),
                    end_date=dt.datetime(2002, 1, 1),
                    name='Salmon',
                    run_name='test',
                    area='4250.6',
                    elevation='843.0',
                    latitude=54.4848,
                    longitude=-123.3659,
                    )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='raven-mohyse',
            datainputs=datainputs)

        assert_response_success(resp)
        tmp_file, _ = urlretrieve(get_output(resp.xml)['hydrograph'])
        ds = xr.open_dataset(tmp_file)
   
        assert ds.variables['q_sim'].shape[0] == 3

