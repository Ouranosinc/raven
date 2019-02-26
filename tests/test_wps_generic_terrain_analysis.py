from pywps import Service
from pywps.tests import assert_response_success
from .common import client_for, TESTDATA, CFG_FILE, get_output

from raven.processes import TerrainAnalysisProcess


class TestGenericTerrainAnalysisProcess:

    def test_simple(self):
        client = client_for(Service(processes=[TerrainAnalysisProcess(), ], cfgfiles=CFG_FILE))
        fields = [
            'raster=file@xlink:href=file://{raster}',
            # 'shape=file@xlink:href=file://{shape}',
            'projected_crs={projected_crs}',
            'band={band}',
            'select_all_touching={touches}',
        ]

        datainputs = ';'.join(fields).format(
            raster=TESTDATA['earthenv_dem_90m'],
            # shape=TESTDATA['donnees_quebec_mrc_poly'],  # TESTDATA['watershed_vector'],
            projected_crs='32198',
            band=1,
            touches=True,
        )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='terrain-analysis', datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)

        assert 'slope' in out
        assert 'aspect' in out
