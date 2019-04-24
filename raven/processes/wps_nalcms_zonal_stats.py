import logging
import json
import tempfile

from pywps import LiteralInput, ComplexInput
from pywps import ComplexOutput
from pywps import Process, FORMATS
from pywps.app.Common import Metadata
from rasterstats import zonal_stats

from raven.utils import archive_sniffer, crs_sniffer, generic_vector_reproject, single_file_check

LOGGER = logging.getLogger("PYWPS")

NALCMS_CATEGORIES = {
    1: 'Forest',
    2: 'Forest',
    3: 'Forest',
    4: 'Forest',
    5: 'Forest',
    6: 'Forest',
    7: 'Shrubs',
    8: 'Shrubs',
    9: 'Grass',
    10: 'Grass',
    11: 'Shrubs',
    12: 'Grass',
    13: 'Grass',
    14: 'Wetland',
    15: 'Crops',
    16: 'Grass',
    17: 'Urban',
    18: 'Water',
    19: 'SnowIce',
}


class NALCMSZonalStatisticsProcess(Process):
    """Given files containing vector data and raster data, perform zonal statistics of the overlapping regions"""

    def __init__(self):
        inputs = [
            ComplexInput('shape', 'Vector Shape',
                         abstract='An ESRI Shapefile, GML, JSON, GeoJSON, or single layer GeoPackage.'
                                  ' The ESRI Shapefile must be zipped and contain the .shp, .shx, and .dbf.'
                                  ' The shape and raster should have a matching CRS.',
                         min_occurs=1, max_occurs=1,
                         supported_formats=[FORMATS.GEOJSON, FORMATS.GML, FORMATS.JSON, FORMATS.SHP]),
            ComplexInput('raster', 'Gridded Land Use raster data set',
                         abstract='The Land Use raster to be queried. Default is the CEC NALCMS (2010)',
                         metadata=[Metadata(
                             'Commission for Environmental Cooperation North American Land Change Monitoring System',
                             'http://www.cec.org/tools-and-resources/map-files/land-cover-2010-landsat-30m'),
                                   Metadata(
                                       "Latifovic, R., Homer, C., Ressl, R., Pouliot, D., Hossain, S.N., Colditz, R.R.,"
                                       " Olthof, I., Giri, C., Victoria, A., 2012. North American land change"
                                       " monitoring system. In: Giri, C., (Ed), Remote Sensing of Land Use and Land"
                                       " Cover: Principles and Applications, CRC-Press, pp. 303-324")],
                         min_occurs=1, max_occurs=1, supported_formats=[FORMATS.GEOTIFF]),
            LiteralInput('return_geojson', 'Return the geometry and statistics as properties in a GeoJSON',
                         data_type='boolean', default='true',
                         min_occurs=1, max_occurs=1),
            LiteralInput('select_all_touching', 'Additionally select boundary pixels that are touched by shape',
                         data_type='boolean', default='false',
                         min_occurs=1, max_occurs=1),
        ]

        outputs = [
            ComplexOutput('statistics', 'DEM properties within the region defined by `shape`.',
                          abstract='Elevation statistics: min, max, mean, median, sum, nodata',
                          supported_formats=[FORMATS.JSON, FORMATS.GEOJSON]),
        ]

        super(NALCMSZonalStatisticsProcess, self).__init__(
            self._handler,
            identifier="nalcms-zonal-stats",
            title="NALCMS Land Use Zonal Statistics",
            version="1.0",
            abstract="Return zonal statistics for the CEC NALCMS based on the boundaries of a vector file.",
            metadata=[],
            inputs=inputs,
            outputs=outputs,
            status_supported=True,
            store_supported=True)

    def _handler(self, request, response):

        # TODO: Deploy the CEC NALCMS 2010 data set on PAVICS GeoServer and access via WCS
        raster_url = request.inputs['raster'][0].file

        shape_url = request.inputs['shape'][0].file
        geojson_out = request.inputs['return_geojson'][0].data
        touches = request.inputs['select_all_touching'][0].data

        vectors = ['.gml', '.shp', '.gpkg', '.geojson', '.json']
        vector_file = single_file_check(archive_sniffer(shape_url, working_dir=self.workdir, extensions=vectors))
        rasters = ['.tiff', '.tif']
        raster_file = single_file_check(archive_sniffer(raster_url, working_dir=self.workdir, extensions=rasters))

        vec_crs, ras_crs = crs_sniffer(vector_file), crs_sniffer(raster_file)

        if vec_crs != ras_crs:
            msg = 'CRS for files {} and {} are not the same. Reprojecting...'.format(vector_file, raster_file)
            LOGGER.warning(msg)

            # Reproject full vector to preserve feature attributes
            projected = tempfile.NamedTemporaryFile(prefix='reprojected_', suffix='.json', delete=False,
                                                    dir=self.workdir).name
            generic_vector_reproject(vector_file, projected, driver='GeoJSON', source_crs=vec_crs, target_crs=ras_crs)
            vector_file = projected

        try:
            stats = zonal_stats(
                vector_file, raster_file, stats=['count', 'nodata'],
                band=1, categorical=True, category_map=NALCMS_CATEGORIES, all_touched=touches,
                geojson_out=geojson_out, raster_out=False)

            if not geojson_out:
                response.outputs['statistics'].data = json.dumps(stats)

            else:
                if len(stats) > 1:
                    feature_collect = {'type': 'FeatureCollection', 'features': stats}
                else:
                    feature_collect = stats

                response.outputs['statistics'].data = json.dumps(feature_collect)

        except Exception as e:
            msg = 'Failed to perform zonal statistics using {} and {}: {}'.format(shape_url, raster_url, e)
            LOGGER.error(msg)
            raise Exception(msg)

        return response
