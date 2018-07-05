from pywps import Process
from pywps import LiteralInput, LiteralOutput
from pywps import ComplexInput, ComplexOutput
from pywps import Format, FORMATS
from pywps.app.Common import Metadata

"""
Dependencies for pysheds not installed with python setup.py install. See requirements.txt. 
numpy, scipy, pandas, geojson, affine, scikit-image, pyproj, rasterio
"""

import logging
LOGGER = logging.getLogger("PYWPS")


class WatershedDelineation(Process):
    """
    Using a DEM and outlet coordinates, return the watershed contour.
    """

    def __init__(self):
        inputs = [
            LiteralInput('latitude', 'Outlet latitude', data_type='float',
                         abstract='Latitudinal coordinate of the watershed outlet.',),
            LiteralInput('longitude', 'Outlet longitude', data_type='float',
                         abstract='Longitudinal coordinate of the watershed outlet.', ),
            LiteralInput('name', 'Watershed name', data_type='string',
                         abstract='Name of the watershed.'),
            ComplexInput('dem', 'Digital Elevation Model',
                         abstract='An URL pointing at the DEM to be used to compute the watershed boundary. Defaults '
                                  'to the HydroSheds DEM.', # TODO: Include details (resolution, version).
                         metadata=[Metadata('HydroSheds Database', 'http://hydrosheds.org'),
                                   Metadata('Lehner, B., Verdin, K., Jarvis, A. (2008): New global hydrography derived from spaceborne elevation data. Eos, Transactions, AGU, 89(10): 93-94.', 'https://doi.org/10.1029/2008EO100001')],
                         min_occurs=0,
                         default='', # TODO: Enter default DEM from PAVICS
                         supported_formats=[FORMATS.GEOTIFF, FORMATS.GML, FORMATS.WCS])
                         ,
        ]
        outputs = [
            ComplexOutput('boundary', 'Watershed boundary',
                          abstract='A polygon defining the watershed boundary.',
                          as_reference=True,
                          supported_formats=[FORMATS.GML]),
        ]

        super(WatershedDelineation, self).__init__(
            self._handler,
            identifier="watershed_delineation",
            title="Watershed delineation algorithm",
            version="1.0",
            abstract="Return the boundary of a watershed computed using a digital elevation model.",
            metadata=[],
            inputs=inputs,
            outputs=outputs,
            status_supported=True,
            store_supported=True)

    @staticmethod
    def _handler(request, response):
        from pysheds.grid import Grid

        # Create pysheds Grid object
        dem_fn = request.inputs['dem'][0].file
        dem = Grid.from_raster(dem_fn, 'dem')

        lat = request.inputs['latitude'][0].data
        lon = request.inputs['longitude'][0].data
        name = request.inputs['name'][0].data

        # Untested
        # dem.catchment(data='dem', x=lon, y=lat, out_name=name, nodata_in=-32768)

        # TODO: Get catchment boundary as a polygon, and convert to adequate output format.
        # TODO: Decide which output formats should be supported.

        # response.outputs['boundary'].data = dem.name

        return response
