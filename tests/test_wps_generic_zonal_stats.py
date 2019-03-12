from pywps import Service
from pywps.tests import assert_response_success
from .common import client_for, TESTDATA, CFG_FILE, get_output

from raven.processes import ZonalStatisticsProcess
import json

from shapely.geometry import shape, MultiPolygon


class TestGenericZonalStatsProcess:

    def test_simple_json(self):
        client = client_for(Service(processes=[ZonalStatisticsProcess(), ], cfgfiles=CFG_FILE))

        fields = [
            'select_all_touching={touches}',
            'return_geojson={return_geojson}',
            'categorical={categorical}',
            'band={band}',
            'shape=file@xlink:href=file://{shape}',
            'raster=file@xlink:href=file://{raster}'
        ]

        datainputs = ';'.join(fields).format(
            touches=True,
            return_geojson=False,
            categorical=False,
            band=1,
            shape=TESTDATA['mrc_subset'],
            raster=TESTDATA['earthenv_dem_90m']
        )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='zonal-stats', datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)

        assert {'statistics'}.issubset([*out])
        assert out['statistics'].endswith('.json')

        with open(out['statistics'], 'r') as f:
            stats = json.load(f)[0]

            assert {'count', 'min', 'max', 'mean', 'median', 'sum', 'nodata'}.issubset([*stats])

    def test_simple_categorized_json(self):
        client = client_for(Service(processes=[ZonalStatisticsProcess(), ], cfgfiles=CFG_FILE))

        fields = [
            'select_all_touching={touches}',
            'return_geojson={return_geojson}',
            'categorical={categorical}',
            'band={band}',
            'shape=file@xlink:href=file://{shape}',
            'raster=file@xlink:href=file://{raster}'
        ]

        datainputs = ';'.join(fields).format(
            touches=True,
            return_geojson=False,
            categorical=True,
            band=1,
            shape=TESTDATA['mrc_subset'],
            raster=TESTDATA['earthenv_dem_90m']
        )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='zonal-stats', datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)

        assert ['statistics'] == [*out]
        assert out['statistics'].endswith('.json')

        with open(out['statistics'], 'r') as f:
            stats = json.load(f)[0]

            assert {'count', 'min', 'max', 'mean', 'median', 'sum', 'nodata'}.issubset([*stats])

            assert set([str(x) for x in range(135, 228)]).issubset([*stats])
            assert '228' not in [*stats]

    def test_simple_geojson(self):
        client = client_for(Service(processes=[ZonalStatisticsProcess(), ], cfgfiles=CFG_FILE))

        fields = [
            'select_all_touching={touches}',
            'return_geojson={return_geojson}',
            'categorical={categorical}',
            'band={band}',
            'shape=file@xlink:href=file://{shape}',
            'raster=file@xlink:href=file://{raster}'
        ]

        datainputs = ';'.join(fields).format(
            touches=True,
            return_geojson=True,
            categorical=False,
            band=1,
            shape=TESTDATA['mrc_subset'],
            raster=TESTDATA['earthenv_dem_90m']
        )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='zonal-stats', datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)

        assert ['statistics'] == [*out]
        assert out['statistics'].endswith('.geojson')

        with open(out['statistics'], 'r') as f:
            features = json.load(f)['features']

            stats = features[0]['properties']
            assert {'count', 'min', 'max', 'mean', 'median', 'sum', 'nodata'}.issubset([*stats])

            geometry = shape(features[0]['geometry'])
            assert isinstance(type(geometry), type(MultiPolygon))

