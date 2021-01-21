import json

import numpy as np
import pytest
from pywps import Service
from pywps.tests import assert_response_success
from ravenpy.utilities.testdata import get_test_data

from raven.processes import HindcastEvaluationProcess

from .common import CFG_FILE, client_for, get_output


class TestForecastEvaluationProcess:
    def test_forecast_evaluation_deterministic(self):
        client = client_for(
            Service(
                processes=[
                    HindcastEvaluationProcess(),
                ],
                cfgfiles=CFG_FILE,
            )
        )

        obs = get_test_data("XSS_forecast_data", "XSS_obs.nc")[0]
        hcst = get_test_data("XSS_forecast_data", "XSS_fcst_det.nc")[0]
        datainputs = (
            f"obs=files@xlink:href=file://{obs};"
            f"hcst=files@xlink:href=file://{hcst};"
            "obs_var=obs;"
            "hcst_var=fcst"
        )
        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="hindcast-evaluation",
            datainputs=datainputs,
        )
        assert_response_success(resp)
        out = get_output(resp.xml)["metrics"]
        m = json.loads(out)
        np.testing.assert_almost_equal(m["mse"], 0.0868818, 4)

    def test_forecast_evaluation_ensemble(self):
        client = client_for(
            Service(
                processes=[
                    HindcastEvaluationProcess(),
                ],
                cfgfiles=CFG_FILE,
            )
        )

        obs = get_test_data("XSS_forecast_data", "XSS_obs.nc")[0]
        hcst = get_test_data("XSS_forecast_data", "XSS_fcst_ens.nc")[0]
        datainputs = (
            f"obs=files@xlink:href=file://{obs};"
            f"hcst=files@xlink:href=file://{hcst};"
            "obs_var=obs;"
            "hcst_var=fcst;"
            "metric=crps_ensemble"
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="hindcast-evaluation",
            datainputs=datainputs,
        )

        assert_response_success(resp)
        out = get_output(resp.xml)["metrics"]
        m = json.loads(out)
        np.testing.assert_almost_equal(m["crps_ensemble"], 0.1791973, 4)
