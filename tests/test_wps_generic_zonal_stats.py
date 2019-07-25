import json

from pywps import Service
from pywps.tests import assert_response_success
from shapely.geometry import shape, MultiPolygon

from raven.processes import ZonalStatisticsProcess
from .common import client_for, TESTDATA, CFG_FILE, get_output


class TestGenericZonalStatsProcess:

    def test_simple(self):
        client = client_for(Service(processes=[ZonalStatisticsProcess(), ], cfgfiles=CFG_FILE))

        fields = [
            'select_all_touching={touches}',
            'categorical={categorical}',
            'band={band}',
            'shape=file@xlink:href=file://{shape}',
            'raster=file@xlink:href=file://{raster}'
        ]

        datainputs = ';'.join(fields).format(
            touches=True,
            categorical=False,
            band=1,
            shape=TESTDATA['mrc_subset'],
            raster=TESTDATA['earthenv_dem_90m']
        )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='zonal-stats', datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)
        feature = json.loads(out['statistics'])['features'][0]

        stats = feature['properties']
        assert {'count', 'min', 'max', 'mean', 'median', 'sum', 'nodata'}.issubset(stats)

        geometry = shape(feature['geometry'])
        assert isinstance(type(geometry), type(MultiPolygon))

    def test_categorized(self):
        client = client_for(Service(processes=[ZonalStatisticsProcess(), ], cfgfiles=CFG_FILE))

        fields = [
            'select_all_touching={touches}',
            'categorical={categorical}',
            'band={band}',
            'shape=file@xlink:href=file://{shape}',
            'raster=file@xlink:href=file://{raster}'
        ]

        datainputs = ';'.join(fields).format(
            touches=True,
            categorical=True,
            band=1,
            shape=TESTDATA['mrc_subset'],
            raster=TESTDATA['earthenv_dem_90m']
        )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='zonal-stats', datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)
        feature = json.loads(out['statistics'])['features'][0]

        stats = feature['properties']
        assert {'count', 'min', 'max', 'mean', 'median', 'sum', 'nodata'}.issubset(stats)

        # Check for accurate pixel counts
        category_count = 0
        for key, val in stats.items():
            try:
                int(key)
                category_count += val
            except ValueError:
                if key in ['count', 'nodata']:
                    continue

        assert (category_count + stats['nodata']) == stats['count']

        geometry = shape(feature['geometry'])
        assert isinstance(type(geometry), type(MultiPolygon))

    def test_geoserver_dem_wcs(self):
        client = client_for(Service(processes=[ZonalStatisticsProcess(), ], cfgfiles=CFG_FILE))
        fields = [
            'select_all_touching={touches}',
            'categorical={categorical}',
            'band={band}',
            'shape=file@xlink:href=file://{shape}', ]

        datainputs = ';'.join(fields).format(
            touches=True,
            categorical=False,
            band=1,
            shape=TESTDATA['mrc_subset'],
        )
        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='zonal-stats', datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)
        feature = json.loads(out['statistics'])['features'][0]

        stats = feature['properties']
        assert {'count', 'min', 'max', 'mean', 'median', 'sum', 'nodata'}.issubset(stats)

        geometry = shape(feature['geometry'])
        assert isinstance(type(geometry), type(MultiPolygon))
