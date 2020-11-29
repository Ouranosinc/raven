import logging
import json
import tempfile

from pywps import LiteralInput, ComplexInput
from pywps import ComplexOutput
from pywps import Process, FORMATS
from pywps.app.Common import Metadata
from rasterstats import zonal_stats
from raven.utils import archive_sniffer, crs_sniffer, single_file_check, generic_vector_reproject
from ravenpy.utilities import gis

LOGGER = logging.getLogger("PYWPS")
SUMMARY_ZONAL_STATS = ['count', 'min', 'max', 'mean', 'median', 'sum', 'nodata']


class ZonalStatisticsProcess(Process):
    """Given files containing vector data and raster data, perform zonal statistics of the overlapping regions"""

    def __init__(self):
        inputs = [
            ComplexInput('shape', 'Vector Shape',
                         abstract='An ESRI Shapefile, GML, JSON, GeoJSON, or single layer GeoPackage.'
                                  ' The ESRI Shapefile must be zipped and contain the .shp, .shx, and .dbf.'
                                  ' The shape and raster should have a matching CRS.',
                         min_occurs=1, max_occurs=1,
                         supported_formats=[FORMATS.GEOJSON, FORMATS.GML, FORMATS.JSON, FORMATS.SHP]),
            ComplexInput('raster', 'Gridded raster data set',
                         abstract='The DEM to be queried. Defaults to the EarthEnv-DEM90 product.',
                         metadata=[Metadata('EarthEnv-DEM90', 'https://www.earthenv.org/DEM'),
                                   Metadata(
                                       'Robinson, Natalie, James Regetz, and Robert P. Guralnick (2014). '
                                       'EarthEnv-DEM90: A Nearly-Global, Void-Free, Multi-Scale Smoothed, 90m Digital '
                                       'Elevation Model from Fused ASTER and SRTM Data. ISPRS Journal of '
                                       'Photogrammetry and Remote Sensing 87: 57–67.',
                                       'https://doi.org/10.1016/j.isprsjprs.2013.11.002')],
                         min_occurs=0, max_occurs=1, supported_formats=[FORMATS.GEOTIFF]),
            LiteralInput('band', 'Raster band',
                         data_type='integer', default=1,
                         abstract='Band of raster examined to perform zonal statistics. Default: 1',
                         min_occurs=1, max_occurs=1),
            LiteralInput('categorical', 'Return distinct pixel categories',
                         data_type='boolean', default='false',
                         min_occurs=1, max_occurs=1),
            LiteralInput('select_all_touching', 'Additionally select boundary pixels that are touched by shape',
                         data_type='boolean', default='false',
                         min_occurs=1, max_occurs=1)
        ]

        outputs = [
            ComplexOutput('statistics', 'DEM properties within the region defined by `shape`.',
                          abstract='Elevation statistics: min, max, mean, median, sum, nodata',
                          supported_formats=[FORMATS.JSON, FORMATS.GEOJSON]),
        ]

        super(ZonalStatisticsProcess, self).__init__(
            self._handler,
            identifier="zonal-stats",
            title="Raster Zonal Statistics",
            version="1.0",
            abstract="Return zonal statistics based on the boundaries of a vector file.",
            metadata=[],
            inputs=inputs,
            outputs=outputs,
            status_supported=True,
            store_supported=True)

    def _handler(self, request, response):

        shape_url = request.inputs['shape'][0].file
        band = request.inputs['band'][0].data
        categorical = request.inputs['categorical'][0].data
        touches = request.inputs['select_all_touching'][0].data

        vectors = ['.gml', '.shp', '.gpkg', '.geojson', '.json']
        vector_file = single_file_check(archive_sniffer(shape_url, working_dir=self.workdir, extensions=vectors))
        rasters = ['.tiff', '.tif']

        if 'raster' in request.inputs:
            raster_url = request.inputs['raster'][0].file
            raster_file = single_file_check(archive_sniffer(raster_url, working_dir=self.workdir, extensions=rasters))
        else:
            bbox = gis.get_bbox(vector_file)
            raster_url = 'public:EarthEnv_DEM90_NorthAmerica'
            raster_bytes = gis.get_raster_wcs(bbox, geographic=True, layer=raster_url)
            raster_file = tempfile.NamedTemporaryFile(prefix='wcs_', suffix='.tiff', delete=False,
                                                      dir=self.workdir).name
            with open(raster_file, 'wb') as f:
                f.write(raster_bytes)

        vec_crs, ras_crs = crs_sniffer(vector_file), crs_sniffer(raster_file)

        if ras_crs != vec_crs:
            msg = 'CRS for files {} and {} are not the same. Reprojecting vector...'.format(vector_file, raster_file)
            LOGGER.warning(msg)

            # Reproject full vector to preserve feature attributes
            projected = tempfile.NamedTemporaryFile(prefix='reprojected_', suffix='.json', delete=False,
                                                    dir=self.workdir).name
            generic_vector_reproject(vector_file, projected, source_crs=vec_crs, target_crs=ras_crs)
            vector_file = projected

        summary_stats = SUMMARY_ZONAL_STATS

        try:
            stats = zonal_stats(
                vector_file, raster_file, stats=summary_stats, band=band, categorical=categorical,
                all_touched=touches, geojson_out=True, raster_out=False)

            feature_collect = {'type': 'FeatureCollection', 'features': stats}
            response.outputs['statistics'].data = json.dumps(feature_collect)

        except Exception as e:
            msg = 'Failed to perform zonal statistics using {} and {}: {}'.format(shape_url, raster_url, e)
            LOGGER.error(msg)
            raise Exception(msg) from e

        return response
