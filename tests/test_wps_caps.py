from common import client_for
from pywps import Service

from raven.processes import processes


def test_wps_caps():
    client = client_for(Service(processes=processes))
    resp = client.get(service="wps", request="getcapabilities", version="1.0.0")
    names = resp.xpath_text(
        "/wps:Capabilities/wps:ProcessOfferings/wps:Process/ows:Identifier"
    )
    sn = set(names.split())

    assert sn == {
        "raster-subset",
        "shape-properties",
        "hydrobasins-select",
        "terrain-analysis",
        "nalcms-zonal-stats",
        "nalcms-zonal-stats-raster",
        "zonal-stats",
    }
