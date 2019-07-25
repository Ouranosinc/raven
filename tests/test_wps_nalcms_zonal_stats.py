import json

from pywps import Service
from pywps.tests import assert_response_success
from shapely.geometry import shape, MultiPolygon

from raven.processes import NALCMSZonalStatisticsProcess
from .common import client_for, TESTDATA, CFG_FILE, get_output


class TestNALCMSZonalStatsProcess:

    def test_simplified_categories(self):
        client = client_for(Service(processes=[NALCMSZonalStatisticsProcess(), ], cfgfiles=CFG_FILE))

        fields = [
            'select_all_touching={touches}',
            'simple_categories={simple_categories}',
            'band={band}',
            'shape=file@xlink:href=file://{shape}',
            'raster=file@xlink:href=file://{raster}'
        ]

        datainputs = ';'.join(fields).format(
            touches=True,
            simple_categories=True,
            band=1,
            shape=TESTDATA['mrc_subset'],
            raster=TESTDATA['cec_nalcms_2010']
        )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='nalcms-zonal-stats', datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)
        feature = json.loads(out['statistics'])['features'][0]

        stats = feature['properties']
        assert {'count', 'nodata'}.issubset(stats)

        category_counts = 0
        for key, val in stats['land-use'].items():
            category_counts += val
        assert category_counts == stats['count']

        geometry = shape(feature['geometry'])
        assert isinstance(type(geometry), type(MultiPolygon))

    def test_true_categories(self):
        client = client_for(Service(processes=[NALCMSZonalStatisticsProcess(), ], cfgfiles=CFG_FILE))

        fields = [
            'select_all_touching={touches}',
            'simple_categories={simple_categories}',
            'band={band}',
            'shape=file@xlink:href=file://{shape}',
            'raster=file@xlink:href=file://{raster}'
        ]

        datainputs = ';'.join(fields).format(
            touches=True,
            simple_categories=False,
            band=1,
            shape=TESTDATA['mrc_subset'],
            raster=TESTDATA['cec_nalcms_2010']
        )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='nalcms-zonal-stats', datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)
        feature = json.loads(out['statistics'])['features'][0]

        stats = feature['properties']
        assert {'count', 'nodata'}.issubset(stats)

        category_counts = 0
        for key, val in stats['land-use'].items():
            category_counts += val
        assert category_counts == stats['count']

        geometry = shape(feature['geometry'])
        assert isinstance(type(geometry), type(MultiPolygon))

    def test_wcs_simplified_categories(self):
        client = client_for(Service(processes=[NALCMSZonalStatisticsProcess(), ], cfgfiles=CFG_FILE))
        fields = [
            'select_all_touching={touches}',
            'simple_categories={simple_categories}',
            'band={band}',
            'shape=file@xlink:href=file://{shape}', ]

        datainputs = ';'.join(fields).format(
            touches=True,
            simple_categories=True,
            band=1,
            shape=TESTDATA['basin_test'],
        )
        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='nalcms-zonal-stats', datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)
        feature = json.loads(out['statistics'])['features'][0]

        stats = feature['properties']
        assert {'count', 'nodata'}.issubset(stats)

        category_counts = 0
        for key, val in stats['land-use'].items():
            category_counts += val
        assert category_counts == stats['count']

        geometry = shape(feature['geometry'])
        assert isinstance(type(geometry), type(MultiPolygon))

    def test_wcs_true_categories(self):
        client = client_for(Service(processes=[NALCMSZonalStatisticsProcess(), ], cfgfiles=CFG_FILE))
        fields = [
            'select_all_touching={touches}',
            'simple_categories={simple_categories}',
            'band={band}',
            'shape=file@xlink:href=file://{shape}', ]

        datainputs = ';'.join(fields).format(
            touches=True,
            simple_categories=False,
            band=1,
            shape=TESTDATA['basin_test'],
        )
        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='nalcms-zonal-stats', datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)
        feature = json.loads(out['statistics'])['features'][0]

        stats = feature['properties']
        assert {'count', 'nodata'}.issubset(stats)

        category_counts = 0
        for key, val in stats['land-use'].items():
            category_counts += val
        assert category_counts == stats['count']

        geometry = shape(feature['geometry'])
        assert isinstance(type(geometry), type(MultiPolygon))
