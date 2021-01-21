from pywps import Service
from pywps.tests import assert_response_success
from ravenpy.utilities.testdata import get_test_data

from raven.processes import RasterSubsetProcess

from .common import CFG_FILE, client_for, get_output


class TestGenericRasterSubsetProcess:
    def test_simple(self):
        client = client_for(
            Service(
                processes=[
                    RasterSubsetProcess(),
                ],
                cfgfiles=CFG_FILE,
            )
        )

        fields = [
            "shape=file@xlink:href=file://{shape}",
            "raster=file@xlink:href=file://{raster}",
            "band={band}",
            "select_all_touching={touches}",
        ]

        datainputs = ";".join(fields).format(
            shape=get_test_data("donneesqc_mrc_poly", "mrc_subset.gml")[0],
            raster=get_test_data(
                "earthenv_dem_90m", "earthenv_dem90_southernQuebec.tiff"
            )[0],
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
