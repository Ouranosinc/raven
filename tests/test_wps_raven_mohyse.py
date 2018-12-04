import pytest
import datetime as dt
import numpy as np

from pywps import Service
from pywps.tests import assert_response_success

from .common import client_for, TESTDATA, CFG_FILE, get_output, urlretrieve
from raven.processes import RavenMOHYSEProcess


class TestRavenMOHYSEProcess:

    def test_simple(self):
        client = client_for(Service(processes=[RavenMOHYSEProcess(), ], cfgfiles=CFG_FILE))

        params = '1.0000, 0.0468, 4.2952, 2.6580, 0.4038, 0.0621, 0.0273, 0.0453, 0.9039, 5.6167'

        datainputs = "ts=files@xlink:href=file://{ts};" \
                     "params={params};" \
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
                    start_date=dt.datetime(1954, 1, 1),
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

        # checking correctness of NSE (full period 1954-2011 would be NSE=0.391107 as template in Wiki)
        assert 'DIAG_NASH_SUTCLIFFE' in tmp_content[0]
        idx_diag = tmp_content[0].split(',').index("DIAG_NASH_SUTCLIFFE")
        diag = np.float(tmp_content[1].split(',')[idx_diag])
        np.testing.assert_almost_equal(diag, 0.3837, 4, err_msg='NSE is not matching expected value')

        # checking correctness of RMSE (full period 1954-2011 would be RMSE=36.7012 as template in wiki)
        assert 'DIAG_RMSE' in tmp_content[0]
        idx_diag = tmp_content[0].split(',').index("DIAG_RMSE")
        diag = np.float(tmp_content[1].split(',')[idx_diag])
        np.testing.assert_almost_equal(diag, 37.1292, 4, err_msg='RMSE is not matching expected value')
