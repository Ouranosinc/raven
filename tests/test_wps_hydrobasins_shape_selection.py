import json

import pytest
from common import CFG_FILE, client_for, get_output
from pywps import Service
from pywps.tests import assert_response_success

from raven.processes import HydroBasinsSelectionProcess
from raven.utils import parse_lonlat


class TestParser:
    def test_parse_lonlat(self):
        lonlats = ["123,345", "122233344.11111 34554554.2", "1111,2222"]
        for ll in lonlats:
            lonlat = parse_lonlat(ll)
            assert len(lonlat) == 2
            assert isinstance(lonlat[0], float)
            assert isinstance(lonlat[1], float)

        with pytest.raises(Exception):
            parse_lonlat("This isn't a number, 333.444")


@pytest.mark.online
class TestShapeSelectionProcess:
    @pytest.mark.very_slow
    @pytest.mark.parametrize("aggregate_upstream,ids", [(True, 551), (False, 1)])
    def test_manicouagan(self, aggregate_upstream, ids):
        client = client_for(
            Service(processes=[HydroBasinsSelectionProcess()], cfgfiles=CFG_FILE)
        )

        fields = [
            "location={location}",
            # 'level={level}',
            # 'lakes={lakes}',
            "aggregate_upstream={aggregate_upstream}",
        ]

        datainputs = ";".join(fields).format(
            location="-68.724444, 50.646667",
            # level="12",
            # lakes=True,
            aggregate_upstream=aggregate_upstream,
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="hydrobasins-select",
            datainputs=datainputs,
        )

        assert_response_success(resp)
        out = get_output(resp.xml)

        assert {"feature", "upstream_ids"}.issubset([*out])
        assert len(json.loads(out["upstream_ids"])) == ids

    @pytest.mark.slow
    def test_lac_saint_jean(self):
        client = client_for(
            Service(processes=[HydroBasinsSelectionProcess()], cfgfiles=CFG_FILE)
        )

        fields = [
            "location={location}",
            # 'level={level}',
            # 'lakes={lakes}',
            "aggregate_upstream={aggregate_upstream}",
        ]

        datainputs_upstream = ";".join(fields).format(
            location="-72.0, 48.5",
            # level="12",
            # lakes=True,
            aggregate_upstream=True,
        )

        datainputs_subbasin = ";".join(fields).format(
            location="-72.0, 48.5",
            # level="12",
            # lakes=True,
            aggregate_upstream=False,
        )

        resp_upstream = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="hydrobasins-select",
            datainputs=datainputs_upstream,
        )
        resp_subbasin = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="hydrobasins-select",
            datainputs=datainputs_subbasin,
        )

        assert_response_success(resp_upstream)
        assert_response_success(resp_subbasin)
        out_upstream = get_output(resp_upstream.xml)
        out_subbasin = get_output(resp_subbasin.xml)

        assert {"feature", "upstream_ids"}.issubset([*out_upstream])
        assert {"feature", "upstream_ids"}.issubset([*out_subbasin])

        assert json.loads(out_subbasin["feature"])["type"] == "FeatureCollection"
        assert json.loads(out_upstream["feature"])["type"] == "Feature"

    @pytest.mark.slow
    def test_smallwood_reservoir(self):
        client = client_for(
            Service(processes=[HydroBasinsSelectionProcess()], cfgfiles=CFG_FILE)
        )

        fields = [
            "location={location}",
            # 'level={level}',
            # 'lakes={lakes}',
            "aggregate_upstream={aggregate_upstream}",
        ]

        datainputs = ";".join(fields).format(
            location="-64.2, 54.1",
            level="12",
            lakes=True,
            aggregate_upstream=True,
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="hydrobasins-select",
            datainputs=datainputs,
        )

        assert_response_success(resp)
        out = get_output(resp.xml)

        assert {"feature", "upstream_ids"}.issubset([*out])

    @pytest.mark.slow
    def test_great_slave_lake(self):
        client = client_for(
            Service(processes=[HydroBasinsSelectionProcess()], cfgfiles=CFG_FILE)
        )

        fields = [
            "location={location}",
            # 'level={level}',
            # 'lakes={lakes}',
            "aggregate_upstream={aggregate_upstream}",
        ]

        datainputs = ";".join(fields).format(
            location="-114.65, 61.35",
            # level="12",
            # lakes=True,
            aggregate_upstream=True,
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="hydrobasins-select",
            datainputs=datainputs,
        )

        assert_response_success(resp)
        out = get_output(resp.xml)

        assert {"feature", "upstream_ids"}.issubset([*out])
