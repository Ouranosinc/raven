import logging
import os

from pywps import LiteralInput, ComplexInput
from pywps import LiteralOutput, ComplexOutput
from pywps import Process, FORMATS
from pywps.app.Common import Metadata
from rasterstats import zonal_stats

from raven.utils import extract_archive

LOGGER = logging.getLogger("PYWPS")


class ZonalStatisticsProcess(Process):
    """Given files containing vector data and raster data, perform zonal statistics of the overlapping regions"""

    identifier = ''
    abstract = ''
    title = ''
    version = ''

    def __init__(self):
        inputs = [
            LiteralInput('select_all_touching', 'Select boundary pixels that are touched by shape',
                         data_type='boolean', default='false'),
            LiteralInput('categorical', 'Return distinct pixel categories',
                         data_type='boolean', default='false'),
            LiteralInput('band', 'Raster band', data_type='integer', default=1,
                         abstract='Band of raster examined to perform zonal statistics. Defaults to 1'),
            ComplexInput('shape', 'Vector Shape',
                         abstract='An URL pointing to either an ESRI Shapefile, GML, GeoJSON, or any other file in a'
                                  ' standard vector format. The ESRI Shapefile must be zipped and contain the .shp,'
                                  ' .shx, and .dbf. The shape CRS definition should also match the DEM CRS.',
                         min_occurs=1, supported_formats=[FORMATS.GEOJSON, FORMATS.GML, FORMATS.JSON, FORMATS.SHP]),
            ComplexInput('raster', 'Gridded Raster Data set',
                         abstract='An URL pointing at the DEM to be queried. Defaults to the USGS HydroSHEDS DEM.',
                         # TODO: Include details (resolution, version).
                         metadata=[Metadata('HydroSheds Database', 'http://hydrosheds.org'),
                                   Metadata(
                                       'Lehner, B., Verdin, K., Jarvis, A. (2008): New global hydrography derived from'
                                       ' spaceborne elevation data. Eos, Transactions, AGU, 89(10): 93-94.',
                                       'https://doi.org/10.1029/2008EO100001')],
                         min_occurs=0, max_occurs=1, supported_formats=[FORMATS.GEOTIFF])]

        outputs = [
            LiteralOutput('count', 'Feature Count', data_type='integer', abstract='Number of features in shape', ),
            LiteralOutput('min', 'Minimum pixel value', data_type='float', abstract='Minimum raster value'),
            LiteralOutput('max', 'Maximum pixel value', data_type='float', abstract='Maximum raster value'),
            LiteralOutput('mean', 'Mean pixel value', data_type='float', abstract='Mean raster value'),
            LiteralOutput('median', 'Median pixel value', data_type='float', abstract='Median raster value'),
            LiteralOutput('sum', 'Sum of pixels', data_type='integer', abstract='Sum of all pixel values'),
            LiteralOutput('nodata', 'Number of null data pixels', data_type='integer',
                          abstract='Number of null data pixels'),
            ComplexOutput('categories', 'Counts of pixels by category',
                          abstract='Counts of pixels by category'),
        ]

        super(ZonalStatisticsProcess, self).__init__(
            self._zonalstats_handler,
            identifier="rasterstats",
            title="Raster Zonal Statistics",
            version="1.0",
            abstract="Return raster zonal statistics based on boundaries of a vector file.",
            metadata=[],
            inputs=inputs,
            outputs=outputs,
            status_supported=True,
            store_supported=True)

    @staticmethod
    def _zonalstats_handler(request, response):

        # dem_fn = request.inputs['raster'][0].file
        # shape_fn = request.inputs['shape'][0].file
        # band = request.inputs['band'][0]
        # touches = request.inputs['select_all_touching'][0]
        # categorical = request.inputs['categorical'][0]

        dem_url = request['raster']
        shape_url = request['shape']
        band = request['band']
        touches = request['select_all_touching']
        categorical = request['categorical']

        archive_types = ['.nc', '.tar', '.zip']
        allowed_types = ['.gml', '.shp', '.geojson', '.json', '.gpkg']

        # TODO: I know this is ugly. Suggestions welcome.
        if any(ext in shape_url for ext in archive_types):
            extracted = extract_archive(shape_url, os.getcwd())  # self.workdir)
            for potential_vector in extracted:
                if any(ext in potential_vector for ext in allowed_types):
                    shape_url = potential_vector

        dem_files = any(dem in extract_archive(dem_url) for dem in ['.tiff', '.tif'])

        print(dem_files)

        try:
            if not categorical:
                stats = zonal_stats(shape_url, dem_url, all_touched=touches, band=band,
                                    stats=['count', 'min', 'max', 'mean', 'median', 'sum', 'nodata'])
                for key in stats.keys():
                    # response.outputs[key].data = stats[key]
                    response[key] = stats[key]
            else:
                stats = zonal_stats(
                    shape_url, dem_url, band=band, categorical=categorical, all_touched=touches, geojson_out=True)
                # response.outputs['categories'].data = stats[1]
                response['categories'] = stats[1]

        except Exception as e:
            msg = 'Failed to perform zonal statistics: {}'.format(e)
            raise Exception(msg)

        return response


if __name__ == "__main__":
    from tests.common import TESTDATA

    fields = ['select_all_touching={select_all_touching}', 'categorical={categorical}', 'band={band}',
              'crs={crs}', 'shape=file@xlink:href=file://{shape}', 'raster=file@xlink:href=file://{raster}']

    inputs = {'select_all_touching': True,
              'categorical': False,
              'band': 1,
              'crs': 4326,
              'shape': TESTDATA['watershed_vector'],
              'raster': TESTDATA['hydrosheds_conditioned']}

    datainputs = ';'.join(fields).format(**inputs)

    response = {'count': None, 'min': None, 'max': None, 'mean': None, 'median': None, 'sum': None, 'nodata': None}

    q = ZonalStatisticsProcess._zonalstats_handler(request=inputs, response=response)

    for k in q.keys():
        print(k, q[k])
