from pywps import Service
import pytest
from pywps.tests import assert_response_success

from .common import client_for
from raven.processes import processes

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
        'sleep',
        'wordcounter']
