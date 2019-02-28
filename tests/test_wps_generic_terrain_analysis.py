from pywps import Service
from pywps.tests import assert_response_success
from .common import client_for, TESTDATA, CFG_FILE, get_output

from raven.processes import TerrainAnalysisProcess
import numpy as np
import rasterio


class TestGenericTerrainAnalysisProcess:

    def test_simple(self):
        client = client_for(Service(processes=[TerrainAnalysisProcess(), ], cfgfiles=CFG_FILE))
        fields = [
            'raster=file@xlink:href=file://{raster}',
            # 'shape=file@xlink:href=file://{shape}',
            'projected_crs={projected_crs}',
            'select_all_touching={touches}',
        ]

        datainputs = ';'.join(fields).format(
            raster=TESTDATA['earthenv_dem_90m'],
            # shape=TESTDATA['watershed_vector'],  # TESTDATA['donnees_quebec_mrc_poly'],
            projected_crs='32198',
            touches=True,
        )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='terrain-analysis', datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)

        assert {'slope', 'aspect'}.issubset([*out])

        with rasterio.open(out['slope'], 'r') as src:
            grid = src.read()
            np.testing.assert_equal(grid.shape, (1, 7279, 5146))
            np.testing.assert_approx_equal(np.max(grid), 68.372086)

        with rasterio.open(out['aspect'], 'r') as src:
            grid = src.read()
            np.testing.assert_equal(grid.shape, (1, 7279, 5146))
            np.testing.assert_approx_equal(np.max(grid), 359.8404)


