import json


import numpy as np
import pytest
from pywps import Service
from pywps.tests import assert_response_success


from raven.processes import ForecastEvaluationProcess
from .common import client_for, TESTDATA, CFG_FILE, get_output



class TestForecastEvaluationProcess:

    def test_forecast_evaluation_deterministic(self):
        client = client_for(Service(processes=[ForecastEvaluationProcess(), ], cfgfiles=CFG_FILE))

        datainputs = "obs=files@xlink:href=file://{obs};" \
                    "fcst=files@xlink:href=file://{fcst};".format(obs=TESTDATA['XSS_obs'],fcst=TESTDATA['XSS_fcst_det'])

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='forecast-evaluation',
            datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)['metrics']
        m = json.loads(out)
        np.testing.assert_almost_equal(m['mse'], 0.0868818, 4)
    
    
    def test_forecast_evaluation_ensemble(self):
        client = client_for(Service(processes=[ForecastEvaluationProcess(), ], cfgfiles=CFG_FILE))

        datainputs = "obs=files@xlink:href=file://{obs};" \
                    "fcst=files@xlink:href=file://{fcst};".format(obs=TESTDATA['XSS_obs'],fcst=TESTDATA['XSS_fcst_ensemble'])

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='forecast-evaluation',
            datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)['metrics']
        m = json.loads(out)
        np.testing.assert_almost_equal(m['crps_ensemble'], 0.1791973, 4)