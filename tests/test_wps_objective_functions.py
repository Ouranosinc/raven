import pytest
import numpy as np
from pywps import Service
from pywps.tests import assert_response_success

from raven.processes import ObjectiveFunctionProcess
from .common import client_for, TESTDATA, CFG_FILE, get_output
from raven.models import Raven
import tempfile
import json


@pytest.fixture(scope="module")
def gr4j():
    rvs = TESTDATA['raven-gr4j-cemaneige-nc-rv']
    ts = TESTDATA['raven-gr4j-cemaneige-nc-ts']

    model = Raven(tempfile.mkdtemp())
    model.configure(rvs)
    model.run([ts, ], )
    return model


class TestObjectiveFunctionProcess:

    def test_all_gr4j(self, gr4j):
        client = client_for(Service(processes=[ObjectiveFunctionProcess(), ], cfgfiles=CFG_FILE))

        kwds = dict(hydrograph=gr4j.outputs['hydrograph'])

        datainputs = "obs=files@xlink:href=file://{hydrograph};"\
                     "sim=files@xlink:href=file://{hydrograph};".format(**kwds)

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='objective-function',
            datainputs=datainputs)

        assert_response_success(resp)
        m = json.loads(get_output(resp.xml)['metrics'])

        np.testing.assert_almost_equal(m['nashsutcliffe'], gr4j.diagnostics['DIAG_NASH_SUTCLIFFE'], 3)

    def test_rmse_gr4j(self, gr4j):
        client = client_for(Service(processes=[ObjectiveFunctionProcess(), ], cfgfiles=CFG_FILE))

        kwds = dict(hydrograph=gr4j.outputs['hydrograph'], name='rmse')

        datainputs = "obs=files@xlink:href=file://{hydrograph};"\
                     "sim=files@xlink:href=file://{hydrograph};"\
                     "name={name}".format(**kwds)

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='objective-function',
            datainputs=datainputs)

        assert_response_success(resp)
        m = json.loads(get_output(resp.xml)['metrics'])

        np.testing.assert_almost_equal(m['rmse'], gr4j.diagnostics['DIAG_RMSE'], 3)
