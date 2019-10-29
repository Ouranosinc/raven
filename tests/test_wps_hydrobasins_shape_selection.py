from pywps import Service
from pywps.tests import assert_response_success
from .common import client_for, CFG_FILE, get_output
import json

from raven.processes import HydroBasinsSelectionProcess


class TestShapeSelectionProcess:

    def test_manicouagan(self):
        client = client_for(Service(processes=[HydroBasinsSelectionProcess(), ], cfgfiles=CFG_FILE))

        fields = [
            'location={location}',
            # 'level={level}',
            # 'lakes={lakes}',
            'aggregate_upstream={aggregate_upstream}',
        ]

        datainputs = ';'.join(fields).format(
            location="-68.724444, 50.646667",
            level="12",
            lakes=True,
            aggregate_upstream=True,
        )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='hydrobasins-select', datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)

        assert {'feature', 'upstream_ids'}.issubset([*out])

    def test_lac_saint_jean(self):
        client = client_for(Service(processes=[HydroBasinsSelectionProcess(), ], cfgfiles=CFG_FILE))

        fields = [
            'location={location}',
            # 'level={level}',
            # 'lakes={lakes}',
            'aggregate_upstream={aggregate_upstream}',
        ]

        datainputs_upstream = ';'.join(fields).format(
            location="-72.0, 48.5",
            level="12",
            lakes=True,
            aggregate_upstream=True,
        )

        datainputs_subbasin = ';'.join(fields).format(
            location="-72.0, 48.5",
            level="12",
            lakes=True,
            aggregate_upstream=False,
        )

        resp_upstream = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='hydrobasins-select',
            datainputs=datainputs_upstream)
        resp_subbasin = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='hydrobasins-select',
            datainputs=datainputs_subbasin)

        assert_response_success(resp_upstream)
        assert_response_success(resp_subbasin)
        out_upstream = get_output(resp_upstream.xml)
        out_subbasin = get_output(resp_subbasin.xml)

        assert {'feature', 'upstream_ids'}.issubset([*out_upstream])
        assert {'feature', 'upstream_ids'}.issubset([*out_subbasin])

        assert json.loads(out_subbasin['feature'])['type'] == 'Feature'
        assert json.loads(out_upstream['feature'])['type'] == 'Feature'

    def test_smallwood_reservoir(self):
        client = client_for(Service(processes=[HydroBasinsSelectionProcess(), ], cfgfiles=CFG_FILE))

        fields = [
            'location={location}',
            # 'level={level}',
            # 'lakes={lakes}',
            'aggregate_upstream={aggregate_upstream}',
        ]

        datainputs = ';'.join(fields).format(
            location="-64.2, 54.1",
            level="12",
            lakes=True,
            aggregate_upstream=True,
        )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='hydrobasins-select', datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)

        assert {'feature', 'upstream_ids'}.issubset([*out])
