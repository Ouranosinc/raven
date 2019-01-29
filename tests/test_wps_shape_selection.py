from pywps import Service
from pywps.tests import assert_response_success
from .common import client_for, TESTDATA, CFG_FILE, get_output

from raven.processes import ShapeSelectionProcess


class TestShapeSelectionProcess:

    def test_simple(self):
        client = client_for(Service(processes=[ShapeSelectionProcess(), ], cfgfiles=CFG_FILE))

        fields = ['collect_upstream={collect_upstream}', 'lonlat_coordinate={lonlat_coordinate}', 'crs={crs}',
                  'shape=file@xlink:href=file://{shape}']
        datainputs = ';'.join(fields).format(
            collect_upstream=True,
            lonlat_coordinate="(-68.724444, 50.646667)",
            crs=4326,
            shape=TESTDATA['hydrobasins_12'],
        )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='shape-selection', datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)

        # TODO: add a couple of explicit tests that properties are computed.

        assert 'properties' in out['properties']
        assert len(out['geometry']) > 1
        assert 'upstream_basins' in out

        # For local testing only; relies on git-lfs stored files
        # assert [7120952471, 7120256602, 7120256491, 7120255681, 7120255891, 7120252371, 7120251081, 7120250931,
        #         7120948591, 7120948041, 7120948581, 7120948491, 7120247371, 7120247861, 7120246521, 7120246421,
        #         7120247921, 7120248001, 7120248771, 7120248632, 7120947481, 7120244611, 7120244451, 7120945251,
        #         7120946701, 7120244801, 7120244601, 7120944781, 7120240071, 7120239952, 7120243661, 7120243651,
        #         7120239011, 7120238831, 7120242601, 7120242591, 7120943201, 7120240661, 7120240522, 7120236361,
        #         7120236292, 7120239431, 7120239441, 7120236681, 7120236501, 7120235601, 7120235491, 7120237881,
        #         7120238051, 7120235481, 7120235351, 7120233771, 7120233621, 7120232791, 7120232781, 7120237521,
        #         7120237751, 7120233901, 7120942181, 7120232281, 7120232191, 7120942261, 7120234431, 7120234671,
        #         7120940851, 7120231001, 7120230852, 7120230101] in out['upstream_basins']
