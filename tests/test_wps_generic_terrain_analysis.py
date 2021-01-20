import json
import pytest

from pywps import Service
from pywps.tests import assert_response_success

from raven.processes import TerrainAnalysisProcess
from .common import client_for, TESTDATA, CFG_FILE, get_output


class TestGenericTerrainAnalysisProcess:

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
            shape=TESTDATA['mrc_subset'],
            projected_crs='6622',
            touches=True,
        )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='terrain-analysis', datainputs=datainputs)

        assert_response_success(resp)
        out = json.loads(get_output(resp.xml)['properties'])

        assert out[0]['elevation'] > 0
        assert out[0]['slope'] > 0
        assert out[0]['aspect'] > 0

    @pytest.mark.online
    @pytest.mark.skip('slow')
    def test_shape_subset_wcs(self):
        client = client_for(Service(processes=[TerrainAnalysisProcess(), ], cfgfiles=CFG_FILE))
        fields = [
            'shape=file@xlink:href=file://{shape}',
            'projected_crs={projected_crs}',
            'select_all_touching={touches}',
        ]

        datainputs = ';'.join(fields).format(
            shape=TESTDATA['mrc_subset'],
            projected_crs='6622',
            touches=True,
        )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='terrain-analysis', datainputs=datainputs)

        assert_response_success(resp)
        out = json.loads(get_output(resp.xml)['properties'])

        assert out[0]['elevation'] > 0
        assert out[0]['slope'] > 0
        assert out[0]['aspect'] > 0


    @pytest.mark.online
    def test_single_polygon(self):
        client = client_for(Service(processes=[TerrainAnalysisProcess(), ], cfgfiles=CFG_FILE))
        fields = [
            'shape=file@xlink:href=file://{shape}',
            'projected_crs={projected_crs}',
            'select_all_touching={touches}',
        ]

        datainputs = ';'.join(fields).format(
            shape=TESTDATA['polygon'],
            projected_crs='6622',
            touches=True,
        )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='terrain-analysis', datainputs=datainputs)

        assert_response_success(resp)
        out = json.loads(get_output(resp.xml)['properties'])

        assert out[0]['elevation'] > 0
        assert out[0]['slope'] > 0
        assert out[0]['aspect'] > 0
