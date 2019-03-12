from pywps import Service
from pywps.tests import assert_response_success
from .common import client_for, TESTDATA, CFG_FILE, get_output

from raven.processes import TerrainAnalysisProcess
import numpy as np
import rasterio
import json


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
        out = json.loads(get_output(resp.xml)['properties'])

        assert out[0]['elevation'] > 0
        assert out[0]['slope'] > 0
        assert out[0]['aspect'] > 0

    def test_shape_subset(self):
        client = client_for(Service(processes=[TerrainAnalysisProcess(), ], cfgfiles=CFG_FILE))
        fields = [
            'raster=file@xlink:href=file://{raster}',
            'shape=file@xlink:href=file://{shape}',
            'projected_crs={projected_crs}',
            'select_all_touching={touches}',
        ]

        datainputs = ';'.join(fields).format(
            raster=TESTDATA['earthenv_dem_90m'],
            shape=TESTDATA['mrc_subset'],  # TESTDATA['donnees_quebec_mrc_poly'],
            projected_crs='32198',
            touches=True,
        )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='terrain-analysis', datainputs=datainputs)

        assert_response_success(resp)
        out = json.loads(get_output(resp.xml)['properties'])

        assert out[0]['elevation'] > 0
        assert out[0]['slope'] > 0
        assert out[0]['aspect'] > 0
