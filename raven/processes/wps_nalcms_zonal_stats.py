import json
import logging
import tempfile

from pywps import ComplexOutput
from pywps import LiteralInput, ComplexInput
from pywps import Process, FORMATS
from pywps.app.Common import Metadata
from rasterstats import zonal_stats

from raven.utilities import gis
from raven.utils import archive_sniffer, crs_sniffer, generic_vector_reproject, single_file_check

LOGGER = logging.getLogger("PYWPS")

categories = {
    'Forest': [1, 2, 3, 4, 5, 6],
    'Shrubs': [7, 8, 11],
    'Grass': [9, 10, 12, 13, 16],
    'Wetland': [14],
    'Crops': [15],
    'Urban': [17],
    'Water': [18],
    'SnowIce': [19]
}

NALCMS_CATEGORIES = {i: cat for (cat, ids) in categories.items() for i in ids}
NALCMS_PROJ4 = '+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs=True'


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
                                 'Latifovic, R., Homer, C., Ressl, R., Pouliot, D., Hossain, S.N., Colditz, R.R.,'
                                 'Olthof, I., Giri, C., Victoria, A., (2012). North American land change '
                                 'monitoring system. In: Giri, C., (Ed), Remote Sensing of Land Use and Land '
                                 'Cover: Principles and Applications, CRC-Press, pp. 303-324')],
                         min_occurs=0, max_occurs=1, supported_formats=[FORMATS.GEOTIFF]),
            LiteralInput('band', 'Raster band',
                         data_type='integer', default=1,
                         abstract='Band of raster examined to perform zonal statistics. Default: 1',
                         min_occurs=1, max_occurs=1),
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

        shape_url = request.inputs['shape'][0].file
        band = request.inputs['band'][0].data
        geojson_out = request.inputs['return_geojson'][0].data
        touches = request.inputs['select_all_touching'][0].data

        vectors = ['.gml', '.shp', '.gpkg', '.geojson', '.json']
        vector_file = single_file_check(archive_sniffer(shape_url, working_dir=self.workdir, extensions=vectors))
        vec_crs = crs_sniffer(vector_file)

        if 'raster' in request.inputs:  # For raster files using the UNFAO Land Cover Classification System (19 types)
            rasters = ['.tiff', '.tif']
            raster_url = request.inputs['raster'][0].file
            raster_file = single_file_check(archive_sniffer(raster_url, working_dir=self.workdir, extensions=rasters))
            ras_crs = crs_sniffer(raster_file)

            if vec_crs != ras_crs:
                msg = 'CRS for files {} and {} are not the same. Reprojecting...'.format(vector_file, raster_file)
                LOGGER.warning(msg)

                # Reproject full vector to preserve feature attributes
                projected = tempfile.NamedTemporaryFile(prefix='reprojected_', suffix='.json', delete=False,
                                                        dir=self.workdir).name
                generic_vector_reproject(vector_file, projected, driver='GeoJSON', source_crs=vec_crs,
                                         target_crs=ras_crs)
                vector_file = projected

        else:  # using the NALCMS data from GeoServer
            projected = tempfile.NamedTemporaryFile(prefix='reprojected_', suffix='.json', delete=False,
                                                    dir=self.workdir).name
            generic_vector_reproject(vector_file, projected, driver='GeoJSON', source_crs=vec_crs,
                                     target_crs=NALCMS_PROJ4)
            vector_file = projected

            bbox = gis.get_bbox(projected)
            raster_url = 'public:EarthEnv_DEM90_NorthAmerica'
            raster_bytes = gis.get_nalcms_wcs(bbox)
            raster_file = tempfile.NamedTemporaryFile(prefix='wcs_', suffix='.tiff', delete=False,
                                                      dir=self.workdir).name
            with open(raster_file, 'wb') as f:
                f.write(raster_bytes)

        try:
            stats = zonal_stats(
                vector_file, raster_file, stats=['count', 'nodata'],
                band=band, categorical=True, category_map=NALCMS_CATEGORIES, all_touched=touches,
                geojson_out=geojson_out, raster_out=False)

            if not geojson_out:
                response.outputs['statistics'].data = json.dumps(stats)

            else:
                feature_collect = {'type': 'FeatureCollection', 'features': stats}
                response.outputs['statistics'].data = json.dumps(feature_collect)

        except Exception as e:
            msg = 'Failed to perform zonal statistics using {} and {}: {}'.format(shape_url, raster_url, e)
            LOGGER.error(msg)
            raise Exception(msg)

        return response
