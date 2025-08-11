import tempfile

import rasterio as rio
from common import CFG_FILE, client_for, get_output
from metalink import download as md
from pywps import Service
from pywps.tests import assert_response_success

from raven.processes import RasterSubsetProcess


class TestGenericRasterSubsetProcess:
    def test_simple(self, yangtze):
        client = client_for(
            Service(processes=[RasterSubsetProcess()], cfgfiles=CFG_FILE)
        )

        fields = [
            "shape=file@xlink:href=file://{shape}",
            "raster=file@xlink:href=file://{raster}",
            "band={band}",
            "select_all_touching={touches}",
        ]

        datainputs = ";".join(fields).format(
            shape=yangtze.fetch("watershed_vector/Basin_test.zip"),
            raster=yangtze.fetch("earthenv_dem_90m/earthenv_dem90_southernQuebec.tiff"),
            band=1,
            touches=True,
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="raster-subset",
            datainputs=datainputs,
        )

        assert_response_success(resp)
        out = get_output(resp.xml)
        assert {"raster"}.issubset([*out])

        raster_dir = md.get(out["raster"], path=tempfile.mkdtemp())
        # Depends on PyMetaLink behaviour?
        assert len(raster_dir) in (1, 5)

        bounds = list()
        for f in raster_dir:
            raster = rio.open(f)
            assert raster.bounds
            bounds.append(raster.bounds)
        assert len(set(bounds)) == len(raster_dir)

    def test_multiple_features_metalink(self, yangtze):
        client = client_for(
            Service(processes=[RasterSubsetProcess()], cfgfiles=CFG_FILE)
        )

        fields = [
            "shape=file@xlink:href=file://{shape}",
            "raster=file@xlink:href=file://{raster}",
            "band={band}",
            "select_all_touching={touches}",
        ]

        datainputs = ";".join(fields).format(
            shape=yangtze.fetch("donneesqc_mrc_poly/mrc_subset.gml"),
            raster=yangtze.fetch("earthenv_dem_90m/earthenv_dem90_southernQuebec.tiff"),
            band=1,
            touches=True,
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="raster-subset",
            datainputs=datainputs,
        )

        assert_response_success(resp)
        out = get_output(resp.xml)
        assert {"raster"}.issubset([*out])

        raster_dir = md.get(out["raster"], path=tempfile.mkdtemp())
        assert len(raster_dir) == 6

        bounds = list()
        for f in raster_dir:
            raster = rio.open(f)
            assert raster.bounds
            bounds.append(raster.bounds)
        assert len(set(bounds)) == len(raster_dir)
