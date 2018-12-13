import json

from pysheds.grid import Grid
from pywps import ComplexInput, ComplexOutput
from pywps import FORMATS
from pywps import LiteralInput
from pywps import Process
from pywps.app.Common import Metadata
from shapely import geometry

# from pywps import LiteralOutput, Format

"""
Dependencies for pysheds not installed with python setup.py install. See requirements.txt.
numpy, scipy, pandas, geojson, affine, scikit-image, pyproj, rasterio
"""

import logging

LOGGER = logging.getLogger("PYWPS")

# N    NE    E    SE    S    SW    W    NW#N    N
# This is hard-coded for now but we'll have to include as a parameter for user provided DEMs.
dirmap = (64, 128, 1, 2, 4, 8, 16, 32)
nodata = -32768


# Think about creating two processes, a HydroShedDelineation process (with options hard-coded)
# and a more general WatershedDelineation.


class WatershedDelineation(Process):
    """
    Using a DEM and outlet coordinates, return the watershed contour.
    """

    def __init__(self):
        inputs = [
            LiteralInput('latitude', 'Outlet latitude', data_type='float',
                         abstract='Latitudinal coordinate of the watershed outlet.', ),
            LiteralInput('longitude', 'Outlet longitude', data_type='float',
                         abstract='Longitudinal coordinate of the watershed outlet.', ),
            LiteralInput('name', 'Watershed name', data_type='string',
                         abstract='Name of the watershed.'),
            ComplexInput('dem', 'Digital Elevation Model',
                         abstract='An URL pointing at the DEM to be used to compute the watershed boundary. Defaults '
                                  'to the HydroSheds DEM.',  # TODO: Include details (resolution, version).
                         metadata=[Metadata('HydroSheds Database', 'http://hydrosheds.org'),
                                   Metadata(
                                       'Lehner, B., Verdin, K., Jarvis, A. (2008): New global hydrography derived from '
                                       'spaceborne elevation data. Eos, Transactions, AGU, 89(10): 93-94.',
                                       'https://doi.org/10.1029/2008EO100001')],
                         min_occurs=0,
                         default='',  # TODO: Enter default DEM from PAVICS
                         supported_formats=[FORMATS.GEOTIFF, FORMATS.GML, FORMATS.WCS]),
            ComplexInput('dir', 'Flow direction grid',
                         abstract='An URL pointing at the flow direction grid to be used to compute the watershed '
                                  'boundary. Defaults to the HydroSheds product. If both the DEM and the flow '
                                  'direction are give, the flow direction supercedes the DEM.',
                         # TODO: Include details (resolution, version).
                         metadata=[Metadata('HydroSheds Database', 'http://hydrosheds.org'),
                                   Metadata(
                                       'Lehner, B., Verdin, K., Jarvis, A. (2008): New global hydrography derived from '
                                       'spaceborne elevation data. Eos, Transactions, AGU, 89(10): 93-94.',
                                       'https://doi.org/10.1029/2008EO100001')],
                         min_occurs=0,
                         default='',  # TODO: Enter default DIR from PAVICS
                         supported_formats=[FORMATS.GEOTIFF, FORMATS.GML, FORMATS.WCS]),

        ]
        outputs = [
            ComplexOutput('boundary', 'Watershed boundary',
                          abstract='A polygon defining the watershed boundary.',
                          as_reference=True,
                          supported_formats=FORMATS.GML),
        ]

        super(WatershedDelineation, self).__init__(
            self._pysheds_handler,
            identifier="watershed_delineation",
            title="Watershed delineation algorithm",
            version="1.0",
            abstract="Return the boundary of a watershed computed using a digital elevation model.",
            metadata=[],
            inputs=inputs,
            outputs=outputs,
            status_supported=True,
            store_supported=True)

    # TODO: David, let me know if you work on this. I'm curious to throw my hat in the ring here.
    @staticmethod
    def _pysheds_handler(request, response):

        # Create pysheds Grid object
        dem_fn = request.inputs['dem'][0].file
        dir_fn = request.inputs['dir'][0].file

        if dem_fn:
            grid = Grid.from_raster(dem_fn, 'dem')
            if not dir_fn:
                grid.flowdir('dem', out_name='dir', nodata=nodata)
            else:
                grid.read_raster(dir_fn, 'dir')
        else:
            grid = Grid.from_raster(dir_fn, 'dir')

        lat = request.inputs['latitude'][0].data
        lon = request.inputs['longitude'][0].data
        # name = request.inputs['name'][0].data

        grid.catchment(data='dir', x=lon, y=lat, dirmap=dirmap, out_name='catch', nodata_in=nodata, xytype='label')

        # This returns multiple polygons, and it's not clear which we should use.
        catch = grid.polygonize(grid.catch.astype('int32'), connectivity=8)
        for p, val in catch:
            if val == 0:  # No idea if this is right or not
                poly = p
                break

        # TODO: Select correct polygon
        # TODO: Decide which output formats should be supported.

        response.outputs['boundary'].data = json.dumps(poly)

        return response

    @staticmethod
    def _saga_handler(request, response):
        # TODO: Refer to https://sourceforge.net/p/saga-gis/wiki/Compiling%20SAGA%20on%20Linux/
        # TODO: Also refer to https://sourceforge.net/p/saga-gis/discussion/790705/thread/4b462702/
        pass


def testing():
    """I understand we need both elevation and drainage direction to delineate the watershed.
    We can compute the drainage direction from the DEM using flowdir.

    Downloaded the au and ca files from hydrosheds at 30sec resolution.
    Issues with affine transformation.

    Flow directions seem to be saved in int16 to limit file size. I assume this is why dirmap takes
    powers of 2 to indicate directions.

    Catchment delineation does not seem to work.
    """

    from matplotlib import pyplot as plt

    dem_fn = '../../tests/testdata/ca_dem_30s/ca_dem_30s/'
    dir_fn = '../../tests/testdata/ca_dir_30s/ca_dir_30s/'

    fv = -32768

    grid = Grid.from_raster(dir_fn, 'dir', nodata=fv)
    grid.read_raster(dem_fn, 'dem', nodata=fv)

    lon, lat = -99.0619, 20.933

    fig, (ax1, ax2) = plt.subplots(1, 2)

    idem = ax1.imshow(grid.view('dem'), extent=grid.extent, cmap='cubehelix', zorder=1, vmin=0)
    plt.colorbar(idem, ax=ax1, label='Elevation (m)')

    idir = ax2.imshow(grid.view('dir'), extent=grid.extent, cmap='viridis', zorder=2, vmin=0)
    boundaries = ([0] + sorted(list(dirmap)))
    plt.colorbar(idir, ax=ax2, boundaries=boundaries, values=sorted(dirmap))

    grid.catchment(data='dir', x=lon, y=lat, dirmap=dirmap, out_name='catch', xytype='label', nodata_in=nodata)
    catch = grid.polygonize(grid.catch.astype('int32'), connectivity=8)
    grid.clip_to('catch')

    for (p, v) in catch:
        poly = geometry.asShape(p)
        ax1.plot(*poly.exterior.xy, color='white')


if __name__ == "__main__":
    testing()
