import pytest
import datetime as dt
import numpy as np

from pywps import Service
from pywps.tests import assert_response_success

from . common import client_for, TESTDATA, CFG_FILE, get_output, urlretrieve
from raven.processes import RavenMultiModelProcess


# @pytest.mark.skip
class TestRavenMultiModelProcess:

    def test_simple(self):
        client = client_for(Service(processes=[RavenMultiModelProcess(), ], cfgfiles=CFG_FILE))

        gr4j = '0.529, -3.396, 407.29, 1.072, 16.9, 0.947'
        hmets = '9.5019, 0.2774, 6.3942, 0.6884, 1.2875, 5.4134, 2.3641, 0.0973, 0.0464, 0.1998, 0.0222, -1.0919, ' \
                '2.6851, 0.3740, 1.0000, 0.4739, 0.0114, 0.0243, 0.0069, 310.7211, 916.1947'

        datainputs = "ts=files@xlink:href=file://{ts};" \
                     "gr4j={gr4j};" \
                     "hmets={hmets};" \
                     "start_date={start_date};" \
                     "end_date={end_date};" \
                     "name={name};" \
                     "area={area};" \
                     "latitude={latitude};" \
                     "longitude={longitude};" \
                     "elevation={elevation};" \
            .format(ts=TESTDATA['raven-gr4j-cemaneige-nc-ts'],
                    gr4j=gr4j,
                    hmets=hmets,
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
            service='WPS', request='Execute', version='1.0.0', identifier='raven-gr4j-cemaneige',
            datainputs=datainputs)

        assert_response_success(resp)

        out = get_output(resp.xml)
        assert 'diagnostics' in out
        tmp_file, _ = urlretrieve(out['diagnostics'])
        tmp_content = open(tmp_file).readlines()

        # checking correctness of NSE (full period 1954-2010 would be NSE=0.511214)
        assert 'DIAG_NASH_SUTCLIFFE' in tmp_content[0]
        idx_diag = tmp_content[0].split(',').index("DIAG_NASH_SUTCLIFFE")
        diag = np.float(tmp_content[1].split(',')[idx_diag])
        np.testing.assert_almost_equal(diag, -0.130318, 4, err_msg='NSE is not matching expected value')

        # checking correctness of RMSE (full period 1954-2010 would be RMSE=32.8827)
        assert 'DIAG_RMSE' in tmp_content[0]
        idx_diag = tmp_content[0].split(',').index("DIAG_RMSE")
        diag = np.float(tmp_content[1].split(',')[idx_diag])
        np.testing.assert_almost_equal(diag, 38.1697, 4, err_msg='RMSE is not matching expected value')
