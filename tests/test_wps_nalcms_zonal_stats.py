from pywps import Service
from pywps.tests import assert_response_success
from .common import client_for, TESTDATA, CFG_FILE, get_output

from raven.processes import NALCMSZonalStatisticsProcess
import json

from shapely.geometry import shape, MultiPolygon


class TestNALCMSZonalStatsProcess:

    def test_simple_categorized_json(self):
        client = client_for(Service(processes=[NALCMSZonalStatisticsProcess(), ], cfgfiles=CFG_FILE))

        fields = [
            'select_all_touching={touches}',
            'return_geojson={return_geojson}',
            'shape=file@xlink:href=file://{shape}',
            'raster=file@xlink:href=file://{raster}'
        ]

        datainputs = ';'.join(fields).format(
            touches=True,
            return_geojson=False,
            shape=TESTDATA['mrc_subset'],
            raster=TESTDATA['cec_nalcms_2010']
        )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='nalcms-zonal-stats', datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)

        stats = json.loads(out['statistics'])[0]
        assert {'count', 'nodata'}.issubset(stats)
        assert {'Forest', 'Shrubs', 'Grass', 'Wetland', 'Crops', 'Urban', 'Water'}.issubset(stats.keys())

    def test_simple_geojson(self):
        client = client_for(Service(processes=[NALCMSZonalStatisticsProcess(), ], cfgfiles=CFG_FILE))

        fields = [
            'select_all_touching={touches}',
            'return_geojson={return_geojson}',
            'shape=file@xlink:href=file://{shape}',
            'raster=file@xlink:href=file://{raster}'
        ]

        datainputs = ';'.join(fields).format(
            touches=True,
            return_geojson=True,
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

        geometry = shape(feature['geometry'])
        assert isinstance(type(geometry), type(MultiPolygon))
