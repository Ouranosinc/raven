import pytest
import datetime as dt
import numpy as np
import xarray as xr

from pywps import Service
from pywps.tests import assert_response_success

from . common import client_for, TESTDATA, CFG_FILE, get_output, urlretrieve
from raven.processes import RavenHMETSProcess


class TestRavenHMETSProcess:

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
            .format(ts=TESTDATA['raven-hmets-nc-ts'],
                    params=params,
                    start_date=dt.datetime(2000, 1, 1),
                    end_date=dt.datetime(2002, 1, 1),
                    init='155,455',
                    name='Salmon',
                    run_name='test-hmets',
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
        np.testing.assert_almost_equal(diag, -2.98165, 4, err_msg='NSE is not matching expected value')

        # checking correctness of RMSE (full period 1954-2011 would be RMSE=28.3759 as template in wiki)
        assert 'DIAG_RMSE' in tmp_content[0]
        idx_diag = tmp_content[0].split(',').index("DIAG_RMSE")
        diag = np.float(tmp_content[1].split(',')[idx_diag])
        np.testing.assert_almost_equal(diag, 71.639, 4, err_msg='RMSE is not matching expected value')

    def test_parallel(self):
        client = client_for(Service(processes=[RavenHMETSProcess(), ], cfgfiles=CFG_FILE))

        params1 = '9.5019, 0.2774, 6.3942, 0.6884, 1.2875, 5.4134, 2.3641, 0.0973, 0.0464, 0.1998, 0.0222, -1.0919, ' \
                  '2.6851, 0.3740, 1.0000, 0.4739, 0.0114, 0.0243, 0.0069, 310.7211, 916.1947'
        params2 = '9.5, 0.3, 6, 0.6, 1.2, 5.4, 2.31, 0.09, 0.04, 0.19, 0.02, -1.09, ' \
                  '2.6, 0.3, 1.0000, 0.4739, 0.0114, 0.0243, 0.0069, 310.7211, 916.1947'
        params3 = '9.5019, 0.2774, 6.3942, 0.6884, 1.2875, 5.4134, 2.3641, 0.0973, 0.0464, 0.1998, 0.0222, -1.0919, ' \
                  '2.0, 0.37, 1.0000, 0.4, 0.01, 0.03, 0.01, 300.0, 920.0'

        # some params in Raven input files are derived from those 21 parameters
        # pdefaults.update({'GR4J_X1_hlf':            pdefaults['GR4J_X1']*1000./2.0})    --> x1 * 1000. / 2.0
        # pdefaults.update({'one_minus_CEMANEIGE_X2': 1.0 - pdefaults['CEMANEIGE_X2']})   --> 1.0 - x6

        datainputs = "ts=files@xlink:href=file://{ts};" \
                     "params={params1};" \
                     "params={params2};" \
                     "params={params3};" \
                     "start_date={start_date};" \
                     "end_date={end_date};" \
                     "name={name};" \
                     "area={area};" \
                     "latitude={latitude};" \
                     "longitude={longitude};" \
                     "elevation={elevation};" \
            .format(ts=TESTDATA['raven-hmets-nc-ts'],
                    params1=params1,
                    params2=params2,
                    params3=params3,
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
            service='WPS', request='Execute', version='1.0.0', identifier='raven-hmets',
            datainputs=datainputs)

        assert_response_success(resp)

        tmp_file, _ = urlretrieve(get_output(resp.xml)['hydrograph'])
        ds = xr.open_dataset(tmp_file)

        assert ds.variables['q_sim'].shape[0] == 3
