from common import client_for
from pywps import Service

from raven.processes import processes

URL = "http://localhost:9099"


def test_describe():
    """Check that owslib can parse the processes' description."""
    from owslib.wps import WebProcessingService

    wps = WebProcessingService(URL, skip_caps=True)
    client = client_for(Service(processes=processes))
    resp = client.get(service="wps", request="getcapabilities", version="1.0.0")
    wps.describeprocess("all", xml=resp.data)
