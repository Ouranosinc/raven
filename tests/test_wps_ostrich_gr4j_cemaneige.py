import pytest
import os
import datetime as dt
import numpy as np

from pywps import Service
from pywps.tests import assert_response_success

from . common import client_for, TESTDATA, CFG_FILE, get_output, urlretrieve
from raven.processes import OstrichGR4JCemaNeigeProcess
import zipfile

# @pytest.mark.skip
class TestOstrichGR4JCemaNeigeProcess:

    def test_simple(self):
        client = client_for(Service(processes=[OstrichGR4JCemaNeigeProcess(), ], cfgfiles=CFG_FILE))

        params = '0.529, -3.396, 407.29, 1.072, 16.9, 0.053'
        lowerBounds = '0.01, -15.0, 10.0, 0.0, 1.0, 0.0'
        upperBounds = '2.5, 10.0, 700.0, 7.0, 30.0, 1.'

        # some params in Raven input files are derived from those 21 parameters
        # pdefaults.update({'GR4J_X1_hlf':            pdefaults['GR4J_X1']*1000./2.0})    --> x1 * 1000. / 2.0
        # pdefaults.update({'one_minus_CEMANEIGE_X2': 1.0 - pdefaults['CEMANEIGE_X2']})   --> 1.0 - x6

        datainputs = "ts=files@xlink:href=file://{ts};" \
                     "algorithm={algorithm};" \
                     "max_iterations={max_iterations};" \
                     "params={params};" \
                     "lowerBounds={lowerBounds};" \
                     "upperBounds={upperBounds};" \
                     "start_date={start_date};" \
                     "duration={duration};" \
                     "name={name};" \
                     "run_name={run_name};" \
                     "area={area};" \
                     "latitude={latitude};" \
                     "longitude={longitude};" \
                     "elevation={elevation};" \
                     "random_seed=0" \
            .format(ts=TESTDATA['ostrich-gr4j-cemaneige-nc-ts'],
                    algorithm='DDS',
                    max_iterations=10,
                    params=params,
                    lowerBounds=lowerBounds,
                    upperBounds=upperBounds,
                    start_date=dt.datetime(1954, 1, 1),
                    duration=208,
                    name='Salmon',
                    run_name='test',
                    area='4250.6',
                    elevation='843.0',
                    latitude=54.4848,
                    longitude=-123.3659,
                    )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='ostrich-gr4j-cemaneige',
            datainputs=datainputs)

        assert_response_success(resp)

        out = get_output(resp.xml)

        assert 'diagnostics' in out
        tmp_file, _ = urlretrieve(out['diagnostics'])
        tmp_content = open(tmp_file).readlines()

        # checking correctness of NSE (full period 1954-2010 with budget of 50 would be NSE=0.5779910)
        assert 'DIAG_NASH_SUTCLIFFE' in tmp_content[0]
        idx_diag = tmp_content[0].split(',').index("DIAG_NASH_SUTCLIFFE")
        diag = np.float(tmp_content[1].split(',')[idx_diag])
        np.testing.assert_almost_equal(diag, 0.486033, 4, err_msg='NSE is not matching expected value')

        # checking correctness of RMSE (full period 1954-2010 with budget of 50 would be RMSE=????)

        assert 'DIAG_RMSE' in tmp_content[0]
        idx_diag = tmp_content[0].split(',').index("DIAG_RMSE")
        diag = np.float(tmp_content[1].split(',')[idx_diag])
        np.testing.assert_almost_equal(diag, 37.1449, 4, err_msg='RMSE is not matching expected value')

        assert "rv_config" in out
        rv_config, _ = urlretrieve(out["rv_config"])
        z = zipfile.ZipFile(rv_config)
        assert len(z.filelist) == 7
