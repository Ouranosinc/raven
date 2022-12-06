import zipfile

from pywps import Service
from pywps.tests import assert_response_success

from raven.processes import RavenProcess

from .common import CFG_FILE, client_for

cf = ["rvi", "rvp", "rvc", "rvh", "rvt"]


class TestRavenProcess:
    def test_gr4j_salmon_nc(self, tmp_path, get_local_testdata):
        client = client_for(Service(processes=[RavenProcess()], cfgfiles=CFG_FILE))

        rvs = get_local_testdata("raven-gr4j-cemaneige/raven-gr4j-salmon.rv?")
        conf = tmp_path / "conf.zip"
        with zipfile.ZipFile(conf, "w") as zf:
            for rv in rvs:
                zf.write(rv, arcname=rv.name)

        ts = get_local_testdata(
            "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
        )
        datainputs = (
            f"ts=files@xlink:href=file://{ts};" f"conf=files@xlink:href=file://{conf}"
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="raven",
            datainputs=datainputs,
        )

        assert_response_success(resp)

    def test_hmets(self, tmp_path, get_local_testdata):
        client = client_for(Service(processes=[RavenProcess()], cfgfiles=CFG_FILE))

        rvs = get_local_testdata("raven-hmets/raven-hmets-salmon.rv?")
        conf = tmp_path / "conf.zip"
        with zipfile.ZipFile(conf, "w") as zf:
            for rv in rvs:
                zf.write(rv, arcname=rv.name)

        ts = get_local_testdata("raven-hmets/Salmon-River-Near-Prince-George_*.rvt")

        datainputs = (
            "ts=files@xlink:href=file://{};"
            "ts=files@xlink:href=file://{};"
            "conf=files@xlink:href=file://{conf}"
        ).format(*ts, conf=conf)

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="raven",
            datainputs=datainputs,
        )
        assert_response_success(resp)

    def test_mohyse(self, tmp_path, get_local_testdata):
        client = client_for(Service(processes=[RavenProcess()], cfgfiles=CFG_FILE))

        rvs = get_local_testdata("raven-mohyse/raven-mohyse-salmon.rv?")
        conf = tmp_path / "conf.zip"
        with zipfile.ZipFile(conf, "w") as zf:
            for rv in rvs:
                zf.write(rv, arcname=rv.name)

        ts = get_local_testdata("raven-mohyse/Salmon-River-Near-Prince-George_*.rvt")

        datainputs = (
            "ts=files@xlink:href=file://{};"
            "ts=files@xlink:href=file://{};"
            "conf=files@xlink:href=file://{conf}"
        ).format(*ts, conf=conf)

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="raven",
            datainputs=datainputs,
        )
        assert_response_success(resp)

    def test_hbv_ec(self, tmp_path, get_local_testdata):
        client = client_for(Service(processes=[RavenProcess()], cfgfiles=CFG_FILE))

        rvs = get_local_testdata("raven-hbvec/raven-hbvec-salmon.rv?")
        conf = tmp_path / "conf.zip"
        with zipfile.ZipFile(conf, "w") as zf:
            for rv in rvs:
                zf.write(rv, arcname=rv.name)

        ts = get_local_testdata("raven-hbvec/Salmon-River-Near-Prince-George_*.rvt")

        datainputs = (
            "ts=files@xlink:href=file://{};"
            "ts=files@xlink:href=file://{};"
            "conf=files@xlink:href=file://{conf}"
        ).format(*ts, conf=conf)

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="raven",
            datainputs=datainputs,
        )
        assert_response_success(resp)
