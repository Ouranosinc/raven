import logging

from pywps import LiteralInput, ComplexInput
from pywps import LiteralOutput
from pywps import Process
from pywps.app.Common import Metadata

LOGGER = logging.getLogger("PYWPS")


class ShapeStatistics(Process):
    """Given a file containing vector data, provide general information and spatial characteristics"""

    def __init__(self):
        inputs = [
            LiteralInput('expand features', 'Examine all features in shapefile?', data_type='boolean'),
            ComplexInput('shape', 'Vector Shape',
                         abstract='An URL pointing to either an ESRI Shapefile, GML, GeoJSON, or any other file in a'
                                  ' standard vector format',
                         min_occurs=1,
                         default=''),
            # ComplexInput('raster', 'Gridded Raster Data set',
            #              abstract='An URL pointing at the DEM to be queried. Defaults to the HydroSheds DEM.',
            #              # TODO: Include details (resolution, version).
            #              metadata=[Metadata('HydroSheds Database', 'http://hydrosheds.org'),
            #                        Metadata(
            #                            'Lehner, B., Verdin, K., Jarvis, A. (2008): New global hydrography derived from '
            #                            'spaceborne elevation data. Eos, Transactions, AGU, 89(10): 93-94.',
            #                            'https://doi.org/10.1029/2008EO100001')],
            #              min_occurs=0,
            #              default='')# TODO: Enter default DEM from PAVICS
        ]

        outputs = [
            LiteralOutput('count', 'Feature Count', data_type='int', abstract='Number of features in shape')
            LiteralOutput('area', 'Area Calculations', data_type='float', abstract='Area of shape in sq. km.'),
            LiteralOutput('centroids', 'Centroid Locations', )
        ]

        super(ShapeStatistics, self).__init__(
            self._shape_stats_handler,
            identifier="shapestats",
            title="Shape Statistics",
            version="1.0",
            abstract="Return general statistics based on a vector file provided.",
            metadata=[],
            inputs=inputs,
            outputs=outputs,
            status_supported=True,
            store_supported=True)

    @ staticmethod
    def _shape_stats_handler(request, response):
        response.append(request)
        return response
