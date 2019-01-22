import pytest
from pywps import Service
from pywps.tests import assert_response_success
from .common import client_for, TESTDATA, CFG_FILE, get_output

from raven.processes import ZonalStatisticsProcess


@pytest.mark.skip
class TestZonalStatsProcess:

    def test_simple(self):
        client = client_for(Service(processes=[ZonalStatisticsProcess(), ], cfgfiles=CFG_FILE))

        fields = [
            'select_all_touching={touches}',
            'return_geometry={return_geometry}'
            'categorical={categorical}',
            'band={band}',
            'shape=file@xlink:href=file://{shape}',
            'raster=file@xlink:href=file://{raster}'
        ]

        datainputs = ';'.join(fields).format(
            touches=True,
            return_geometry=False,
            categorical=False,
            band=1,
            shape=TESTDATA['watershed_vector'],
            raster=TESTDATA['earthenv_dem_90m']
        )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='raster-stats', datainputs=datainputs)
        print(resp.response[0])
        assert_response_success(resp)

        out = get_output(resp.xml)

        assert 'properties' in out

        # There is a bug in pywps 4.0 such that as_reference=False is not respected when ComplexOutput.file is set.
        # So the following won't work out of the box.

        # out = out['properties']

        # assert 'count' in out
        # assert 'min' in out
        # assert 'max' in out
        # assert 'mean' in out
        # assert 'median' in out
        # assert 'sum' in out
        # assert 'nodata' in out
        # assert 'categories' not in out
