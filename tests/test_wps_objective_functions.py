import json
import tempfile

import numpy as np
import pytest
from pywps import Service
from pywps.tests import assert_response_success
from ravenpy.models import Raven
from ravenpy.utilities.testdata import get_test_data

from raven.processes import GraphObjectiveFunctionFitProcess, ObjectiveFunctionProcess

from .common import CFG_FILE, client_for, get_output


@pytest.fixture(scope="module")
def gr4j():
    rvs = get_test_data("raven-gr4j-cemaneige", "raven-gr4j-salmon.rv?")
    ts = get_test_data(
        "raven-gr4j-cemaneige", "Salmon-River-Near-Prince-George_meteo_daily.nc"
    )[0]

    model = Raven(tempfile.mkdtemp())
    model.configure(rvs)
    model([ts])
    return model


class TestObjectiveFunctionProcess:
    def test_all_gr4j(self, gr4j):
        client = client_for(
            Service(
                processes=[
                    ObjectiveFunctionProcess(),
                ],
                cfgfiles=CFG_FILE,
            )
        )

        kwds = dict(hydrograph=gr4j.outputs["hydrograph"])

        datainputs = (
            "obs=files@xlink:href=file://{hydrograph};"
            "sim=files@xlink:href=file://{hydrograph};".format(**kwds)
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="objective-function",
            datainputs=datainputs,
        )

        assert_response_success(resp)
        out = get_output(resp.xml)["metrics"]
        m = json.loads(out)

        np.testing.assert_almost_equal(
            m["nashsutcliffe"], gr4j.diagnostics["DIAG_NASH_SUTCLIFFE"], 4
        )

    def test_rmse_gr4j(self, gr4j):
        client = client_for(
            Service(
                processes=[
                    ObjectiveFunctionProcess(),
                ],
                cfgfiles=CFG_FILE,
            )
        )

        kwds = dict(hydrograph=gr4j.outputs["hydrograph"], name="rmse")

        datainputs = (
            "obs=files@xlink:href=file://{hydrograph};"
            "sim=files@xlink:href=file://{hydrograph};"
            "name={name}".format(**kwds)
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="objective-function",
            datainputs=datainputs,
        )

        assert_response_success(resp)
        out = get_output(resp.xml)["metrics"]
        m = json.loads(out)

        np.testing.assert_almost_equal(m["rmse"], gr4j.diagnostics["DIAG_RMSE"], 4)

    def test_wps_graph_objective_function_fit(self, gr4j):
        client = client_for(
            Service(
                processes=[
                    GraphObjectiveFunctionFitProcess(),
                ],
                cfgfiles=CFG_FILE,
            )
        )

        datainputs = "sims=files@xlink:href=file://{sims};".format(
            sims=gr4j.outputs["hydrograph"]
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="graph_objective_function_fit",
            datainputs=datainputs,
        )

        assert_response_success(resp)
