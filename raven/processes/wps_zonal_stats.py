import logging
import os
import json
from pywps import LiteralInput, ComplexInput
from pywps import LiteralOutput, ComplexOutput
from pywps import Process, FORMATS
from pywps.app.Common import Metadata
from rasterstats import zonal_stats

from raven.utils import extract_archive

LOGGER = logging.getLogger("PYWPS")


class ZonalStatisticsProcess(Process):
    """Given files containing vector data and raster data, perform zonal statistics of the overlapping regions"""

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
            ComplexOutput('properties', 'DEM properties within the region defined by `shape`.',
                          abstract='Elevation statistics: min, max, mean, median, sum, nodata',
                          supported_formats=(FORMATS.JSON,),
                          as_reference=False),
        ]
        #            LiteralOutput('count', 'Feature Count', data_type='integer', abstract='Number of features in
        #            shape', ),
        #            LiteralOutput('min', 'Minimum pixel value', data_type='float', abstract='Minimum raster value'),
        #            LiteralOutput('max', 'Maximum pixel value', data_type='float', abstract='Maximum raster value'),
        #            LiteralOutput('mean', 'Mean pixel value', data_type='float', abstract='Mean raster value'),
        #            LiteralOutput('median', 'Median pixel value', data_type='float', abstract='Median raster value'),
        #            LiteralOutput('sum', 'Sum of pixels', data_type='integer', abstract='Sum of all pixel values'),
        #            LiteralOutput('nodata', 'Number of null data pixels', data_type='integer',
        #                          abstract='Number of null data pixels'),

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

        # dem_fn = request.inputs['raster'][0].file
        # shape_fn = request.inputs['shape'][0].file
        # band = request.inputs['band'][0]
        # touches = request.inputs['select_all_touching'][0]
        # categorical = request.inputs['categorical'][0]

        raster_url = request.inputs['raster'][0].file
        shape_url = request.inputs['shape'][0].file
        band = request.inputs['band'][0].data
        touches = request.inputs['select_all_touching'][0].data
        categorical = request.inputs['categorical'][0].data

        archive_types = ['.nc', '.tar', '.zip']
        allowed_vector = ['.gml', '.shp', '.geojson', '.json', '.gpkg']
        allowed_raster = ['.tiff', '.tif']

        vector_file = []
        raster_file = []

        # TODO: I know this is ugly. Suggestions welcome.
        if any(ext in shape_url for ext in archive_types):
            extracted = extract_archive(shape_url, self.workdir)
            for potential_vector in extracted:
                if any(ext in potential_vector for ext in allowed_vector):
                    vector_file = potential_vector

        if any(dem in str(raster_url) for dem in archive_types):
            extracted = extract_archive(raster_url, os.getcwd())  # self.workdir)
            for potential_raster in extracted:
                if any(ext in potential_raster for ext in allowed_raster):
                    raster_file = potential_raster

        try:
            if not categorical:
                stats = zonal_stats(vector_file, raster_file, all_touched=touches, band=band,
                                    stats=['count', 'min', 'max', 'mean', 'median', 'sum', 'nodata'])

            else:
                stats = zonal_stats(
                    shape_url, raster_url, band=band, categorical=categorical, all_touched=touches, geojson_out=True)

            # Using the PyWPS 4.0 release this will output garbage. Should be fixed for the next version.
            # response.outputs['properties'].data = json.dumps(stats)

            # For the time being, this should do, but will return a reference, not the actual content.
            fn = os.path.join(self.workdir, 'prop.json')
            with open(fn, 'w') as f:
                json.dump(stats, f)
            response.outputs['properties'].file = fn

        except Exception as e:
            msg = 'Failed to perform zonal statistics: {}'.format(e)
            raise Exception(msg)

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
#     response = {}  # {'count': None, 'min': None, 'max': None, 'mean': None, 'median': None, 'sum': None, 'nodata': None}
#
#     q = ZonalStatisticsProcess._zonalstats_handler(request=inputs, response=response)
#
#     for k in q.keys():
#         print(k, q[k])
