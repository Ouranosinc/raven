from pywps import Service
from pywps.tests import assert_response_success
from ravenpy.utilities.testdata import get_test_data

from raven.processes import RavenProcess

from .common import CFG_FILE, client_for

cf = ["rvi", "rvp", "rvc", "rvh", "rvt"]


class TestRavenProcess:
    def test_gr4j_salmon_nc(self):
        client = client_for(
            Service(
                processes=[
                    RavenProcess(),
                ],
                cfgfiles=CFG_FILE,
            )
        )

        rvs = get_test_data("raven-gr4j-cemaneige", "raven-gr4j-salmon.rv?")
        ts = get_test_data(
            "raven-gr4j-cemaneige", "Salmon-River-Near-Prince-George_meteo_daily.nc"
        )[0]
        config = {f.suffix[1:]: f for f in rvs}

        datainputs = (
            "ts=files@xlink:href=file://{ts};"
            + ";".join(["conf=files@xlink:href=file://{%s}" % key for key in cf])
        ).format(ts=ts, **config)

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="raven",
            datainputs=datainputs,
        )

        assert_response_success(resp)

    def test_hmets(self):
        client = client_for(
            Service(
                processes=[
                    RavenProcess(),
                ],
                cfgfiles=CFG_FILE,
            )
        )

        rvs = get_test_data("raven-hmets", "raven-hmets-salmon.rv?")
        ts = get_test_data("raven-hmets", "Salmon-River-Near-Prince-George_*.rvt")
        config = {f.suffix[1:]: f for f in rvs}

        datainputs = (
            "ts=files@xlink:href=file://{};"
            "ts=files@xlink:href=file://{};"
            + ";".join(["conf=files@xlink:href=file://{%s}" % key for key in cf])
        ).format(*ts, **config)

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="raven",
            datainputs=datainputs,
        )
        assert_response_success(resp)

    def test_mohyse(self):
        client = client_for(
            Service(
                processes=[
                    RavenProcess(),
                ],
                cfgfiles=CFG_FILE,
            )
        )

        rvs = get_test_data("raven-mohyse", "raven-mohyse-salmon.rv?")
        ts = get_test_data("raven-mohyse", "Salmon-River-Near-Prince-George_*.rvt")
        config = {f.suffix[1:]: f for f in rvs}

        datainputs = (
            "ts=files@xlink:href=file://{};"
            "ts=files@xlink:href=file://{};"
            + ";".join(["conf=files@xlink:href=file://{%s}" % key for key in cf])
        ).format(*ts, **config)

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="raven",
            datainputs=datainputs,
        )
        assert_response_success(resp)

    def test_hbv_ec(self):
        client = client_for(
            Service(
                processes=[
                    RavenProcess(),
                ],
                cfgfiles=CFG_FILE,
            )
        )

        rvs = get_test_data("raven-hbv-ec", "raven-hbv-ec-salmon.rv?")
        ts = get_test_data("raven-hbv-ec", "Salmon-River-Near-Prince-George_*.rvt")
        config = {f.suffix[1:]: f for f in rvs}

        datainputs = (
            "ts=files@xlink:href=file://{};"
            "ts=files@xlink:href=file://{};"
            + ";".join(["conf=files@xlink:href=file://{%s}" % key for key in cf])
        ).format(*ts, **config)

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="raven",
            datainputs=datainputs,
        )
        assert_response_success(resp)
