from urllib.request import urlretrieve

import numpy as np
import xarray as xr
from pywps import Service
from pywps.tests import assert_response_success
from ravenpy.utilities.testdata import get_local_testdata

from raven.processes import ForecastFloodRiskProcess

from .common import CFG_FILE, client_for, get_output


class TestForecastEvaluationProcess:
    def test_forecast_floodrisk_deterministic(self):
        client = client_for(
            Service(processes=[ForecastFloodRiskProcess()], cfgfiles=CFG_FILE)
        )

        datainputs = (
            "fcst=files@xlink:href=file://{fcst};"
            "name=fcst;"
            "flood_level={flood_level};".format(
                fcst=get_local_testdata("flood_risk/XSS_fcst_det.nc"),
                flood_level=0.5,
            )
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="forecast-floodrisk",
            datainputs=datainputs,
        )

        assert_response_success(resp)
        out = get_output(resp.xml)
        floodrisk, _ = urlretrieve(out["flood_risk"])
        ds = xr.open_dataset(floodrisk)
        np.testing.assert_almost_equal(ds["fcst"], [1, 0, 0], 2)

    def test_forecast_floodrisk_ensemble(self):
        client = client_for(
            Service(
                processes=[
                    ForecastFloodRiskProcess(),
                ],
                cfgfiles=CFG_FILE,
            )
        )

        datainputs = (
            "fcst=files@xlink:href=file://{fcst};"
            "name=fcst;"
            "flood_level={flood_level};".format(
                fcst=get_local_testdata("flood_risk/XSS_fcst_ens.nc"), flood_level=0.5
            )
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="forecast-floodrisk",
            datainputs=datainputs,
        )

        assert_response_success(resp)
        out = get_output(resp.xml)
        floodrisk, _ = urlretrieve(out["flood_risk"])
        ds = xr.open_dataset(floodrisk)
        np.testing.assert_almost_equal(ds["fcst"], [0.2, 0.2, 0.8], 2)
