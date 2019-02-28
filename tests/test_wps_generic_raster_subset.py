from pywps import Service
from pywps.tests import assert_response_success
from .common import client_for, TESTDATA, CFG_FILE, get_output

from raven.processes import RasterSubsetProcess


class TestGenericRasterSubsetProcess:

    def test_simple(self):
        client = client_for(Service(processes=[RasterSubsetProcess(), ], cfgfiles=CFG_FILE))

        fields = [
            'shape=file@xlink:href=file://{shape}',
            'raster=file@xlink:href=file://{raster}',
            'band={band}',
            'select_all_touching={touches}',
        ]

        datainputs = ';'.join(fields).format(
            shape=TESTDATA['mrc_subset'],
            raster=TESTDATA['earthenv_dem_90m'],
            band=1,
            touches=True,
        )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='raster-subset', datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)

        assert {'raster'}.issubset([*out])
