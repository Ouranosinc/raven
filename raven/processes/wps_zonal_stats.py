import logging
import os
import json
from pywps import LiteralInput, ComplexInput
from pywps import ComplexOutput
from pywps import Process, FORMATS
from pywps.app.Common import Metadata
from rasterstats import zonal_stats

from raven.utils import archive_sniffer

LOGGER = logging.getLogger("PYWPS")


class ZonalStatisticsProcess(Process):
    """Given files containing vector data and raster data, perform zonal statistics of the overlapping regions"""

    def __init__(self):
        inputs = [
            LiteralInput('select_all_touching', 'Select boundary pixels that are touched by shape',
                         data_type='boolean', default='false'),
            LiteralInput('return_geometry', 'Return the vertices of the geometry in the JSON file ',
                         data_type='boolean', default='false'),
            LiteralInput('categorical', 'Return distinct pixel categories',
                         data_type='boolean', default='false'),
            LiteralInput('band', 'Raster band', data_type='integer', default=1,
                         abstract='Band of raster examined to perform zonal statistics. Defaults to 1'),
            ComplexInput('shape', 'Vector Shape',
                         abstract='An URL pointing to either an ESRI Shapefile, GML, GeoJSON, or any other file in a'
                                  ' standard vector format. The ESRI Shapefile must be zipped and contain the .shp,'
                                  ' .shx, and .dbf. The shape and raster should have a matching CRS.',
                         min_occurs=1, supported_formats=[FORMATS.GEOJSON, FORMATS.GML, FORMATS.JSON, FORMATS.SHP]),
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
            ComplexOutput('properties', 'DEM properties within the region defined by `shape`.',
                          abstract='Elevation statistics: min, max, mean, median, sum, nodata',
                          supported_formats=(FORMATS.JSON, FORMATS.GEOJSON),
                          as_reference=False),
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
        touches = request.inputs['select_all_touching'][0].data
        geojson_out = request.inputs['return_geometry'][0].data
        categorical = request.inputs['categorical'][0].data

        types = ['.tar', '.zip', '.7z']
        vectors = ['.gml', '.shp', '.geojson', '.json']  # '.gpkg' requires more handling
        rasters = ['.tiff', '.tif']

        vector_file = archive_sniffer(shape_url, working_dir=self.workdir, archive_types=types, extensions=vectors)
        raster_file = archive_sniffer(raster_url, working_dir=self.workdir, archive_types=types, extensions=rasters)

        if not raster_file:
            raster_file = raster_url

        try:
            stats = zonal_stats(
                vector_file, raster_file, stats=['count', 'min', 'max', 'mean', 'median', 'sum', 'nodata'],
                band=band, categorical=categorical, all_touched=touches, geojson_out=geojson_out)

            # Using the PyWPS 4.0 release this will output garbage. Should be fixed for the next version.
            # response.outputs['properties'].data = json.dumps(stats)

            # For the time being, this should do, but will return a reference, not the actual content.
            fn = os.path.join(self.workdir, 'prop.json')
            with open(fn, 'w') as f:
                json.dump(stats, f)
            response.outputs['properties'].file = fn

        except Exception as e:
            msg = 'Failed to perform zonal statistics: {}'.format(e)
            LOGGER.error(msg)

        return response


# if __name__ == "__main__":
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
#
#     response = {}
#
#     q = ZonalStatisticsProcess._zonalstats_handler(request=inputs, response=response)
#
#     for k in q.keys():
#         print(k, q[k])
