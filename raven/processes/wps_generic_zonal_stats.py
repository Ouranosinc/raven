import logging
import json
import tempfile

import numpy as np
import rasterio as rio

from affine import Affine
from fiona.crs import from_epsg
from pywps import LiteralInput, ComplexInput
from pywps import ComplexOutput
from pywps import Process, FORMATS
from pywps.app.Common import Metadata
from rasterstats import zonal_stats


from raven.utils import archive_sniffer, crs_sniffer

LOGGER = logging.getLogger("PYWPS")


class ZonalStatisticsProcess(Process):
    """Given files containing vector data and raster data, perform zonal statistics of the overlapping regions"""

    def __init__(self):
        inputs = [
            LiteralInput('select_all_touching', 'Select boundary pixels that are touched by shape',
                         data_type='boolean', default='false'),
            LiteralInput('return_geojson', 'Return the vertices of the geometry in a GeoJSON',
                         data_type='boolean', default='true'),
            LiteralInput('return_raster', 'Return the masked raster in a GeoTIFF',
                         data_type='boolean', default='false'),
            LiteralInput('categorical', 'Return distinct pixel categories',
                         data_type='boolean', default='true'),
            LiteralInput('band', 'Raster band', data_type='integer', default=1,
                         abstract='Band of raster examined to perform zonal statistics. Default: 1'),
            LiteralInput('crs', 'Coordinate Reference System of shape (EPSG) and lat/lon coordinates',
                         data_type='integer',
                         default=4326),
            ComplexInput('shape', 'Vector Shape',
                         abstract='An URL pointing to either an ESRI Shapefile, GML, GeoJSON, or any other file in a'
                                  ' standard vector format. The ESRI Shapefile must be zipped and contain the .shp,'
                                  ' .shx, and .dbf. The shape and raster should have a matching CRS.',
                         min_occurs=1, max_occurs=1,
                         supported_formats=[FORMATS.GEOJSON, FORMATS.GML, FORMATS.JSON, FORMATS.SHP]),
            ComplexInput('raster', 'Gridded Raster Data set',
                         abstract='An URL pointing at the DEM to be queried. Defaults to the USGS HydroSHEDS DEM.',
                         # TODO: Include details (resolution, version).
                         metadata=[Metadata('HydroSheds Database', 'http://hydrosheds.org'),
                                   Metadata(
                                       'Lehner, B., Verdin, K., Jarvis, A. (2008): New global hydrography derived from'
                                       ' spaceborne elevation data. Eos, Transactions, AGU, 89(10): 93-94.',
                                       'https://doi.org/10.1029/2008EO100001')],
                         min_occurs=1, max_occurs=1, supported_formats=[FORMATS.GEOTIFF])]

        outputs = [
            ComplexOutput('statistics', 'DEM properties within the region defined by `shape`.',
                          abstract='Elevation statistics: min, max, mean, median, sum, nodata',
                          supported_formats=[FORMATS.JSON, FORMATS.GEOJSON]),
            ComplexOutput('raster', 'DEM subset of `shape` region in GeoTIFF format.',
                          abstract='Elevation statistics: min, max, mean, median, sum, nodata',
                          supported_formats=[FORMATS.JSON, FORMATS.GEOJSON]),
        ]

        super(ZonalStatisticsProcess, self).__init__(
            self._handler,
            identifier="raster-stats",
            title="Raster Zonal Statistics",
            version="1.0",
            abstract="Return raster zonal statistics based on boundaries of a vector file.",
            metadata=[],
            inputs=inputs,
            outputs=outputs,
            status_supported=True,
            store_supported=True)

    def _handler(self, request, response):

        shape_url = request.inputs['shape'][0].file
        raster_url = request.inputs['raster'][0].file
        band = request.inputs['band'][0].data
        crs = from_epsg(request.inputs['crs'][0].data)
        touches = request.inputs['select_all_touching'][0].data
        geojson_out = request.inputs['return_geojson'][0].data
        raster_out = request.inputs['return_raster'][0].data
        categorical = request.inputs['categorical'][0].data

        vectors = ['.gml', '.shp', '.geojson', '.json']  # '.gpkg' requires more handling
        rasters = ['.tiff', '.tif']

        vector_file = archive_sniffer(shape_url, working_dir=self.workdir, extensions=vectors)
        raster_file = archive_sniffer(raster_url, working_dir=self.workdir, extensions=rasters)

        raster_subset = ''

        crs = crs_sniffer(vector_file, crs)
        try:
            stats = zonal_stats(
                vector_file, raster_file, stats=['count', 'min', 'max', 'mean', 'median', 'sum', 'nodata'],
                band=band, categorical=categorical, all_touched=touches, geojson_out=geojson_out, raster_out=raster_out)

            # Using the PyWPS 4.0 release this will output garbage. Should be fixed for the next version.
            # response.outputs['properties'].data = json.dumps(stats)

            shape_stats = tempfile.NamedTemporaryFile(suffix='.json', delete=False)

            if not raster_out:
                try:
                    with open(shape_stats.name, 'w') as f:
                        json.dump(stats[0], f)
                except Exception as e:
                    msg = 'Failed to write geojson to {}: {}'.format(shape_stats, e)
                    LOGGER.error(msg)

            else:
                if geojson_out:
                    raster_location = stats[0]['properties']
                else:
                    raster_location = stats[0]

                try:
                    raster = raster_location['mini_raster_array']
                    grid_properties = raster_location['mini_raster_affine'][0:6]
                    nodata = raster_location['mini_raster_nodata']

                    aff = Affine(*grid_properties)

                    raster_subset = tempfile.NamedTemporaryFile(suffix='.tiff', delete=False)
                    LOGGER.warning('Writing raster data to {}'.format(raster_subset.name))

                    masked_array = np.ma.masked_values(raster, nodata)
                    data_type = rio.dtypes.get_minimum_dtype(np.max(masked_array))

                    # Make unsigned data types become signed ones!
                    if data_type.startswith('u'):
                        data_type = str(data_type)[1:]
                    normal_array = np.asarray(masked_array, dtype=data_type)

                    with rio.open(raster_subset.name, 'w', driver='GTiff', count=1, height=raster.shape[0],
                                  width=raster.shape[1], dtype=data_type, transform=aff, crs=crs, nodata=nodata) as f:
                        f.write(normal_array, 1)

                    del raster_location['mini_raster_array']
                    del raster_location['mini_raster_affine']
                    del raster_location['mini_raster_nodata']
                    with open(shape_stats.name, 'w') as f:
                        json.dump(stats[0], f)

                except Exception as e:
                    msg = 'Failed to write raster to {}: {}'.format(raster_subset, e)
                    LOGGER.error(msg)

            response.outputs['statistics'].data = shape_stats.name
            response.outputs['raster'].data = raster_subset.name

        except Exception as e:
            msg = 'Failed to perform zonal statistics using {} and {}: {}'.format(shape_url, raster_url, e)
            LOGGER.error(msg)

        return response

# def testing():
#
#     from tests.common import TESTDATA
#
#     fields = ['select_all_touching={select_all_touching}', 'categorical={categorical}', 'band={band}',
#               'crs={crs}', 'shape=file@xlink:href=file://{shape}', 'raster=file@xlink:href=file://{raster}']
#
#     inputs = {'select_all_touching': True,
#               'categorical': True,
#               'band': 1,
#               'crs': 4326,
#               'shape': TESTDATA['watershed_vector'],
#               'raster': TESTDATA['hydrosheds_conditioned']}
#
#     datainputs = ';'.join(fields).format(**inputs)
#     outputs = {}
#
#     q = ZonalStatisticsProcess._handler(request=inputs, response=outputs)
#
#     for k in q.keys():
#         print(k, q[k])
#
#
# if __name__ == "__main__":
#     testing()
