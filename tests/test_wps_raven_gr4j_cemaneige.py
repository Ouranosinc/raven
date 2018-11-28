import pytest
import datetime as dt
import numpy as np

from pywps import Service
from pywps.tests import assert_response_success

from . common import client_for, TESTDATA, CFG_FILE, get_output, urlretrieve
from raven.processes import RavenGR4JCemaNeigeProcess


class TestRavenGR4JCemaNeigeProcess:

    def test_simple(self):
        client = client_for(Service(processes=[RavenGR4JCemaNeigeProcess(), ], cfgfiles=CFG_FILE))

        params = '0.696, 0.7, 19.7, 2.09, 123.3, 0.25'

        # some params in Raven input files are derived from those 21 parameters
        # pdefaults.update({'GR4J_X1_hlf':            pdefaults['GR4J_X1']*1000./2.0})    --> x1 * 1000. / 2.0
        # pdefaults.update({'one_minus_CEMANEIGE_X2': 1.0 - pdefaults['CEMANEIGE_X2']})   --> 1.0 - x6

        datainputs = "ts=files@xlink:href=file://{ts};" \
                     "params={params};" \
                     "start_date={start_date};" \
                     "end_date={end_date};" \
                     "init={init};" \
                     "name={name};" \
                     "area={area};" \
                     "latitude={latitude};" \
                     "longitude={longitude};" \
                     "elevation={elevation};" \
            .format(ts=TESTDATA['raven-gr4j-cemaneige-nc-ts'],
                    params=params,
                    start_date=dt.datetime(2000, 1, 1),
                    end_date=dt.datetime(2002, 1, 1),
                    init='4,5',
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

        # checking correctness of NSE (full period 1954-2011 would be NSE=-0.738921)
        assert 'DIAG_NASH_SUTCLIFFE' in tmp_content[0]
        idx_diag = tmp_content[0].split(',').index("DIAG_NASH_SUTCLIFFE")
        diag = np.float(tmp_content[1].split(',')[idx_diag])
        np.testing.assert_almost_equal(diag, -0.0905, 4, err_msg='NSE is not matching expected value')

        # checking correctness of RMSE (full period 1954-2011 would be RMSE=62.0224)
        assert 'DIAG_RMSE' in tmp_content[0]
        idx_diag = tmp_content[0].split(',').index("DIAG_RMSE")
        diag = np.float(tmp_content[1].split(',')[idx_diag])
        np.testing.assert_almost_equal(diag, 37.5124, 4, err_msg='RMSE is not matching expected value')
