import pytest
from pywps import Service
from pywps.tests import assert_response_success
from ravenpy.utilities.testdata import get_test_data

from raven.processes import GR4JCemaNeigeProcess

from .common import CFG_FILE, client_for


@pytest.mark.skip
class TestGR4JCemaNeigeProcess:
    def test_simple(self):
        client = client_for(
            Service(
                processes=[
                    GR4JCemaNeigeProcess(),
                ],
                cfgfiles=CFG_FILE,
            )
        )
        pr = get_test_data("gr4j_cemaneige", "pr.nc")[0]
        tas = get_test_data("gr4j_cemaneige", "tas.nc")[0]
        evap = get_test_data("gr4j_cemaneige", "evap.nc")[0]

        datainputs = (
            f"pr=files@xlink:href=file://{pr};"
            f"tas=files@xlink:href=file://{tas};"
            f"evap=files@xlink:href=file://{evap};"
        )
        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="gr4j-cemaneige",
            datainputs=datainputs,
        )

        assert_response_success(resp)
