import pytest
from pywps import Service

from raven.processes import processes
from .common import client_for


def test_wps_caps():
    client = client_for(Service(processes=processes))
    resp = client.get(service='wps', request='getcapabilities', version='1.0.0')
    names = resp.xpath_text('/wps:Capabilities'
                            '/wps:ProcessOfferings'
                            '/wps:Process'
                            '/ows:Identifier')
    sn = set(names.split())
    assert sn == {'gr4j-cemaneige',
                  'raster-stats',
                  'raven',
                  'raven-gr4j-cemaneige',
                  'raven-mohyse',
                  'raven-hmets',
                  'raven-hbv-ec',
                  'shape-area',
                  'shape-selection',
                  'objective-function',
                  # 'regionalisation',
                  }
