from pywps import Service
from pywps.tests import assert_response_success
from .common import client_for, TESTDATA, CFG_FILE, get_output

from raven.processes import RasterSubsetProcess


class TestGenericRasterSubsetProcess:

    def test_simple(self):
        client = client_for(Service(processes=[RasterSubsetProcess(), ], cfgfiles=CFG_FILE))

        fields = [
            'select_all_touching={touches}',
            'band={band}',
            'crs={crs}',
            'shape=file@xlink:href=file://{shape}',
            'raster=file@xlink:href=file://{raster}'
        ]

        datainputs = ';'.join(fields).format(
            touches=True,
            band=1,
            crs=4326,
            # shape=TESTDATA['donnees_quebec_mrc_poly'],
            # shape=TESTDATA['watershed_vector'],
            shape=TESTDATA['broken_data'],
            raster=TESTDATA['earthenv_dem_90m']
        )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='raster-subset', datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)

        print(out)

        assert 'raster' in out
