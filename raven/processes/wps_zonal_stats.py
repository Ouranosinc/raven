import logging

from pywps import LiteralInput, ComplexInput
from pywps import LiteralOutput
from pywps import Process, FORMATS
from pywps.app.Common import Metadata
from rasterstats import zonal_stats

LOGGER = logging.getLogger("PYWPS")


class ZonalStatistics(Process):
    """Given a file containing vector data, provide general information and spatial characteristics"""

    def __init__(self):
        inputs = [
            LiteralInput('select_all_touching', 'Select boundary pixels that are touched by shape',
                         data_type='boolean', default='false'),
            LiteralInput('categorical', 'Return distinct pixel categories',
                         data_type='boolean', default='false'),
            LiteralInput('band', 'Raster band', data_type='int', default=1,
                         abstract='Band of raster examined to perform zonal statistics. Defaults to 1'),
            ComplexInput('shape', 'Vector Shape',
                         abstract='An URL pointing to either an ESRI Shapefile, GML, GeoJSON, or any other file in a'
                                  ' standard vector format. The ESRI Shapefile must be zipped and contain the .shp,'
                                  ' .shx, and .dbf. The shape CRS definition should also match the DEM CRS.',
                         min_occurs=1, supported_formats=[FORMATS.GEOJSON, FORMATS.GML, FORMATS.JSON, FORMATS.SHP]),
            # supported_formats=[
            #            # kml
            #            {mimeType: 'text/xml',
            #            encoding: 'windows-1250',
            #            schema: 'http://schemas.opengis.net/kml/2.2.0/ogckml22.xsd'}
            #            ])
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
            LiteralOutput('count', 'Feature Count', data_type='int', abstract='Number of features in shape', ),
            LiteralOutput('min', 'Minimum pixel value', data_type='float', abstract='Minimum raster value'),
            LiteralOutput('max', 'Maximum pixel value', data_type='float', abstract='Maximum raster value'),
            LiteralOutput('mean', 'Mean pixel value', data_type='float', abstract='Mean raster value'),
            LiteralOutput('median', 'Median pixel value', data_type='float', abstract='Median raster value'),
            LiteralOutput('sum', 'Sum of pixels', data_type='int', abstract='Sum of all pixel values'),
            LiteralOutput('nodata', 'Number of null data pixels', data_type='int',
                          abstract='Number of null data pixels'),
            LiteralOutput('categories', 'Counts of pixels by category', data_type='dict',
                          abstract='Counts of pixels by category'),
        ]

        super(ZonalStatistics, self).__init__(
            self._rasterstats_handler,
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
    def _rasterstats_handler(request, response):

        dem_fn = request.inputs['raster'][0].file
        shape_fn = request.inputs['shape'][0].file
        band = request.inputs['band'][0]
        touches = request.inputs['select_all_touching'][0]
        categorical = request.inputs['categorical'][0]

        # TODO: Figure out whether lists are accepted outputs
        try:
            if not categorical:
                stats = zonal_stats(shape_fn, dem_fn, all_touched=touches, band=band,
                                    stats=['count', 'min', 'max', 'mean', 'median', 'sum', 'nodata'])
                for key in stats.keys():
                    response.outputs[key].data = stats[key]

            else:
                stats = zonal_stats(
                    shape_fn, dem_fn, band=band, categorical=categorical, all_touched=touches, geojson_out=True)
                response.outputs['categories'].data = stats[1]

        except Exception as e:
            msg = 'Failed to perform zonal statistics: {}'.format(e)
            raise Exception(msg)

        return response
