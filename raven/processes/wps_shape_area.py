import logging
from functools import partial
from raven.utils import address_append

import fiona
import shapely.ops as ops
from fiona.crs import from_epsg, from_string
from pyproj import Proj, transform
from pywps import LiteralInput, ComplexInput
from pywps import LiteralOutput
from pywps import Process, FORMATS
from shapely.geometry import shape

LOGGER = logging.getLogger("PYWPS")

LAEA = '+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
WORLDMOLL = '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
ALBERS = '+proj=aea +lat_1=0 +lat_2=0 +lon_0=0 +ellps=WGS84 +x_0=0 +y_0=0 +units=m +no_defs'


class ShapeAreaProcess(Process):
    """Given a file containing vector data, provide general information and spatial characteristics"""

    identifier = ''
    abstract = ''
    title = ''
    version = ''

    def __init__(self):
        inputs = [
            LiteralInput('use_all_features', 'Examine all features in shapefile', data_type='boolean', default='false'),
            LiteralInput('crs', 'EPSG Coordinate Reference System code', data_type='integer', default='4326'),
            ComplexInput('shape', 'Vector Shape',
                         abstract='An URL pointing to either an ESRI Shapefile, GML, JSON, GeoJSON.'
                                  ' The ESRI Shapefile must be zipped and contain the .shp, .shx, and .dbf.',
                         min_occurs=1, supported_formats=[FORMATS.GML, FORMATS.GEOJSON, FORMATS.SHP, FORMATS.JSON])]
        # supported_formats=[
        #            # kml
        #            {mimeType: 'text/xml',
        #            encoding: 'windows-1250',
        #            schema: 'http://schemas.opengis.net/kml/2.2.0/ogckml22.xsd'}
        #            ])

        outputs = [
            LiteralOutput('area', 'Area Calculations', data_type='float', abstract='Area of shape in m^2'),
            LiteralOutput('centroids', 'Centroid Locations', data_type='float',
                          abstract="Geographic locations of feature centroids"),
            LiteralOutput('schemas', 'Feature schemas',
                          abstract='Geographic representations and descriptions of shape features')
        ]

        super(ShapeAreaProcess, self).__init__(
            self._shapearea_handler,
            identifier="shape_area",
            title="Shape Area",
            version="1.0",
            abstract="Return shape area in square Kilometres based on line boundaries of a polygonal vector file.",
            metadata=[],
            inputs=inputs,
            outputs=outputs,
            status_supported=True,
            store_supported=True)

    @staticmethod
    def _shapearea_handler(request, response):

        shape_fn = request.inputs['shape'][0].file
        shape_url = address_append(shape_fn)
        shape_crs = request.inputs['crs'][0]
        features = request.inputs['use_all_features'][0]
        schemas = []
        centroids = []
        areas = []

        for i, layername in enumerate(fiona.listlayers(shape_url)):
            if features:
                with fiona.open(shape_url, 'r', crs=from_epsg(shape_crs), layer=i) as src:
                    geom = shape(src['geometry'])

                    schemas.append(src.schema)
                    centroids.append(geom.centroid)
                    transformed = ops.transform(
                        partial(
                            transform,
                            Proj(src.crs),
                            Proj(from_string(LAEA))),
                        geom)
                    areas.append(transformed.area)

            else:
                with fiona.open(shape_url, 'r', crs=from_epsg(shape_crs)) as src:
                    geom = shape(src['geometry'])

                    schemas.append(src.schema)
                    centroids.append(geom.centroid)
                    transformed = ops.transform(
                        partial(
                            transform,
                            Proj(src.crs),
                            Proj(from_string(LAEA))),
                        geom)
                    areas.append(transformed.area)

            src.close()

        response.outputs['area'].data = areas
        response.outputs['centroids'].data = centroids
        response.outputs['schemas'].data = schemas

        return response
