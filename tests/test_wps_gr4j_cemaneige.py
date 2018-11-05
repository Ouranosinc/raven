import pytest

from pywps import Service
from pywps.tests import assert_response_success

from . common import client_for, TESTDATA, CFG_FILE
from raven.processes import GR4JCemaNeigeProcess


class TestGR4JCemaNeigeProcess:

    def test_simple(self):
        client = client_for(Service(processes=[GR4JCemaNeigeProcess(), ], cfgfiles=CFG_FILE))
        kwds = TESTDATA['gr4j-cemaneige']
        datainputs = "pr=files@xlink:href=file://{pr};tas=files@xlink:href=file://{tas};evap=files@xlink:href=file://{evap};".format(**kwds)
        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='gr4j_cemaneige',
            datainputs=datainputs)
        print(datainputs)
        print(resp.response)
        assert_response_success(resp)
