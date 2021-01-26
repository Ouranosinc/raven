import pytest
from pywps import Service
from pywps.tests import assert_response_success
from ravenpy.utilities.testdata import get_local_testdata

from raven.processes import GraphFcstUncertaintyProcess

from .common import CFG_FILE, client_for, get_output


class TestForecastGraphProcess:
    def test_forecast_graph(self):
        client = client_for(
            Service(
                processes=[
                    GraphFcstUncertaintyProcess(),
                ],
                cfgfiles=CFG_FILE,
            )
        )

        datainputs = (
            "fcst=files@xlink:href=file://{fcst};"
            "fcst_var={fcst_var};".format(
                fcst=get_local_testdata("flood_risk", "XSS_fcst_ens.nc"),
                fcst_var="fcst",
            )
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="graph_forecast_uncertainty",
            datainputs=datainputs,
        )

        assert_response_success(resp)
        out = get_output(resp.xml)
        assert out["graph_forecasts"].endswith(".png")

    def test_hindcast_graph(self):
        client = client_for(
            Service(
                processes=[
                    GraphFcstUncertaintyProcess(),
                ],
                cfgfiles=CFG_FILE,
            )
        )

        datainputs = (
            "fcst=files@xlink:href=file://{fcst};"
            "qobs=files@xlink:href=file://{qobs};"
            "fcst_var={fcst_var};"
            "qobs_var={qobs_var};".format(
                fcst=get_local_testdata("flood_risk/XSS_fcst_ens.nc"),
                fcst_var="fcst",
                qobs=get_local_testdata("XSS_forecast_data/XSS_obs.nc"),
                qobs_var="obs",
            )
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="graph_forecast_uncertainty",
            datainputs=datainputs,
        )

        assert_response_success(resp)
        out = get_output(resp.xml)
        assert out["graph_forecasts"].endswith(".png")
