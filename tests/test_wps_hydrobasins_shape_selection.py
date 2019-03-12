from pywps import Service
from pywps.tests import assert_response_success
from .common import client_for, CFG_FILE, get_output

from raven.processes import ShapeSelectionProcess


class TestShapeSelectionProcess:

    def test_manicouagan(self):
        client = client_for(Service(processes=[ShapeSelectionProcess(), ], cfgfiles=CFG_FILE))

        fields = [
            'lonlat_coordinate={lonlat_coordinate}',
            'level={level}',
            'lakes={lakes}',
            'collect_upstream={collect_upstream}',
        ]

        datainputs = ';'.join(fields).format(
            lonlat_coordinate="(-68.724444, 50.646667)",
            level="12",
            lakes=True,
            collect_upstream=True,
        )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='shape-selection', datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)

        assert {'geojson', 'upstream_basins'}.issubset([*out])

    def test_lac_saint_jean(self):
        client = client_for(Service(processes=[ShapeSelectionProcess(), ], cfgfiles=CFG_FILE))

        fields = [
            'lonlat_coordinate={lonlat_coordinate}',
            'level={level}',
            'lakes={lakes}',
            'collect_upstream={collect_upstream}',
        ]

        datainputs_upstream = ';'.join(fields).format(
            lonlat_coordinate="(-72.0, 48.5)",
            level="12",
            lakes=True,
            collect_upstream=True,
        )

        datainputs_subbasin = ';'.join(fields).format(
            lonlat_coordinate="(-72.0, 48.5)",
            level="12",
            lakes=True,
            collect_upstream=False,
        )

        resp_upstream = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='shape-selection',
            datainputs=datainputs_upstream)
        resp_subbasin = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='shape-selection',
            datainputs=datainputs_subbasin)

        assert_response_success(resp_upstream)
        assert_response_success(resp_subbasin)
        out_upstream = get_output(resp_upstream.xml)
        out_subbasin = get_output(resp_subbasin.xml)

        assert {'geojson', 'upstream_basins'}.issubset([*out_upstream])
        assert {'geojson', 'upstream_basins'}.issubset([*out_subbasin])

    def test_smallwood_reservoir(self):
        client = client_for(Service(processes=[ShapeSelectionProcess(), ], cfgfiles=CFG_FILE))

        fields = [
            'lonlat_coordinate={lonlat_coordinate}',
            'level={level}',
            'lakes={lakes}',
            'collect_upstream={collect_upstream}',
        ]

        datainputs = ';'.join(fields).format(
            lonlat_coordinate="(-64.2, 54.1)",
            level="12",
            lakes=True,
            collect_upstream=True,
        )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='shape-selection', datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)

        assert {'geojson', 'upstream_basins'}.issubset([*out])
