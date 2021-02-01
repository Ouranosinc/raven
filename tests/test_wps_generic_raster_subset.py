import tempfile

from metalink import download as md
import rasterio as rio
from ravenpy.utilities.testdata import get_local_testdata

from pywps import Service
from pywps.tests import assert_response_success
from raven.processes import RasterSubsetProcess

from .common import CFG_FILE, client_for, get_output


class TestGenericRasterSubsetProcess:
    def test_simple(self):
        client = client_for(
            Service(
                processes=[
                    RasterSubsetProcess()
                ],
                cfgfiles=CFG_FILE
            )
        )

        fields = [
            "shape=file@xlink:href=file://{shape}",
            "raster=file@xlink:href=file://{raster}",
            "band={band}",
            "select_all_touching={touches}",
        ]

        datainputs = ";".join(fields).format(
            shape=get_local_testdata("donneesqc_mrc_poly/mrc_subset.gml"),
            raster=get_local_testdata(
                "earthenv_dem_90m/earthenv_dem90_southernQuebec.tiff"
            ),
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

        bounds = list()
        for f in raster_dir:
            raster = rio.open(f)
            assert raster.bounds
            bounds.append(raster.bounds)
        assert len(set(bounds)) == len(bounds)
