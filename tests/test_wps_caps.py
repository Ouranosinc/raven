import pytest
from pywps import Service

from raven.processes import processes
from .common import client_for


# from pywps.tests import assert_response_success

@pytest.mark.skip
def test_wps_caps():
    client = client_for(Service(processes=processes))
    resp = client.get(service='wps', request='getcapabilities', version='1.0.0')
    names = resp.xpath_text('/wps:Capabilities'
                            '/wps:ProcessOfferings'
                            '/wps:Process'
                            '/ows:Identifier')
    assert sorted(names.split()) == [
        'inout',
        'gr4j-cemaneige',
        'raven',
        'raven-gr4j-cemaneige',
        'raven-hmets',
        'sleep',
        'wordcounter',]
