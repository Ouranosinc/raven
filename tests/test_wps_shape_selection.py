from pywps import Service
from pywps.tests import assert_response_success
from .common import client_for, TESTDATA, CFG_FILE, get_output

from raven.processes import ShapeSelectionProcess


class TestShapeSelectionProcess:

    def test_simple(self):
        client = client_for(Service(processes=[ShapeSelectionProcess(), ], cfgfiles=CFG_FILE))

        fields = ['lonlat_coordinate={lonlat_coordinate}', 'crs={crs}', 'shape=file@xlink:href=file://{file}']
        datainputs = ';'.join(fields).format(
            lonlat_coordinate="(-68.724444, 50.646667)",
            crs=4326,
            file=TESTDATA['hydrobasins_12'],
        )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='shape-selection', datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)

        assert 'feature' in out
        # TODO: add a couple of explicit tests that properties are computed.
