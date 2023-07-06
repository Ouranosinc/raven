import json

import pytest
from pywps import Service
from pywps.tests import assert_response_success

from raven.processes import TerrainAnalysisProcess

from .common import CFG_FILE, client_for, get_output


class TestGenericTerrainAnalysisProcess:
    def test_shape_subset(self, get_local_testdata):
        client = client_for(
            Service(processes=[TerrainAnalysisProcess()], cfgfiles=CFG_FILE)
        )
        fields = [
            "raster=file@xlink:href=file://{raster}",
            "shape=file@xlink:href=file://{shape}",
            "projected_crs={projected_crs}",
            "select_all_touching={touches}",
        ]

        datainputs = ";".join(fields).format(
            raster=get_local_testdata(
                "earthenv_dem_90m/earthenv_dem90_southernQuebec.tiff"
            ),
            shape=get_local_testdata("donneesqc_mrc_poly/mrc_subset.gml"),
            projected_crs="6622",
            touches=True,
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="terrain-analysis",
            datainputs=datainputs,
        )

        assert_response_success(resp)
        out = json.loads(get_output(resp.xml)["properties"])

        assert out[0]["elevation"] > 0
        assert out[0]["slope"] > 0
        assert out[0]["aspect"] > 0

    @pytest.mark.slow
    def test_shape_subset_wcs(self, get_local_testdata):
        client = client_for(
            Service(processes=[TerrainAnalysisProcess()], cfgfiles=CFG_FILE)
        )
        fields = [
            "shape=file@xlink:href=file://{shape}",
            "projected_crs={projected_crs}",
            "select_all_touching={touches}",
        ]

        datainputs = ";".join(fields).format(
            shape=get_local_testdata("donneesqc_mrc_poly/mrc_subset.gml"),
            projected_crs="6622",
            touches=True,
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="terrain-analysis",
            datainputs=datainputs,
        )

        assert_response_success(resp)
        out = json.loads(get_output(resp.xml)["properties"])

        assert out[0]["elevation"] > 0
        assert out[0]["slope"] > 0
        assert out[0]["aspect"] > 0

    @pytest.mark.online
    def test_single_polygon(self, get_local_testdata):
        client = client_for(
            Service(processes=[TerrainAnalysisProcess()], cfgfiles=CFG_FILE)
        )
        fields = [
            "shape=file@xlink:href=file://{shape}",
            "projected_crs={projected_crs}",
            "select_all_touching={touches}",
        ]

        datainputs = ";".join(fields).format(
            shape=get_local_testdata("polygons/Basin_10.zip"),
            projected_crs="6622",
            touches=True,
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="terrain-analysis",
            datainputs=datainputs,
        )

        assert_response_success(resp)
        out = json.loads(get_output(resp.xml)["properties"])

        assert out[0]["elevation"] > 0
        assert out[0]["slope"] > 0
        assert out[0]["aspect"] > 0
