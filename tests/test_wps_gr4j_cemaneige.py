import pytest
from pywps import Service
from pywps.tests import assert_response_success

from raven.processes import GR4JCemaNeigeProcess
from .common import client_for, TESTDATA, CFG_FILE


@pytest.mark.skip
class TestGR4JCemaNeigeProcess:

    def test_simple(self):
        client = client_for(Service(processes=[GR4JCemaNeigeProcess(), ], cfgfiles=CFG_FILE))
        kwds = TESTDATA['gr4j-cemaneige']
        datainputs = "pr=files@xlink:href=file://{pr};"\
                     "tas=files@xlink:href=file://{tas};"\
                     "evap=files@xlink:href=file://{evap};".format(**kwds)
        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='gr4j-cemaneige',
            datainputs=datainputs)

        assert_response_success(resp)
