from pywps import Service
from pywps.tests import assert_response_success
from .common import client_for, CFG_FILE, get_output

from raven.processes import ShapeSelectionProcess


class TestShapeSelectionProcess:

    def test_simple(self):
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

        # TODO: add a couple of explicit tests that properties are computed.

        assert {'geojson', 'upstream_basins'}.issubset([*out])
