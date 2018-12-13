import logging
from functools import partial
from raven.utils import extract_archive

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
ALBERS_NAM = '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs'


class ShapeAreaProcess(Process):
    """Given a file containing vector data, provide general information and spatial characteristics"""

    def __init__(self):
        inputs = [
            ComplexInput('shape', 'Vector Shape',
                         abstract='An URL pointing to either an ESRI Shapefile, GML, JSON, GeoJSON.'
                                  ' The ESRI Shapefile must be zipped and contain the .shp, .shx, and .dbf.',
                         min_occurs=1,
                         supported_formats=[FORMATS.GML, FORMATS.GEOJSON, FORMATS.SHP, FORMATS.JSON]),
            LiteralInput('crs', 'Coordinate Reference System of shape (EPSG)',
                         data_type='integer',
                         default=4326),
            LiteralInput('projected_crs',
                         'Coordinate Reference System for area calculation (EPSG; Default:32198)',
                         data_type='integer',
                         default=32198)
        ]

        outputs = [
            LiteralOutput('properties', 'Feature schemas',
                          abstract='Geographic representations and descriptions of shape features'),
            LiteralOutput('area', 'Area calculation in square metres', data_type='float',
                          abstract='Area of shape in m^2'),
            LiteralOutput('centroid_lon', 'Centroid longitude',
                          data_type='float',
                          abstract="Geographic location of feature's centroid"),
            LiteralOutput('centroid_lat', 'Centroid latitude',
                          data_type='float',
                          abstract="Geographic location of feature's centroid")
        ]

        super(ShapeAreaProcess, self).__init__(
            self._handler,
            identifier="shape-area",
            title="Shape Area",
            version="1.0",
            abstract="Return shape area in square metres based on line boundaries of a polygonal vector file.",
            metadata=[],
            inputs=inputs,
            outputs=outputs,
            status_supported=True,
            store_supported=True)

    def _handler(self, request, response):

        shape_url = request.inputs['shape'][0].file
        shape_crs = request.inputs['crs'][0].data
        projected_crs = request.inputs['projected_crs'][0].data

        archive_types = ['.nc', '.tar', '.zip']
        allowed_types = ['.gml', '.shp', '.geojson', '.json', '.gpkg']

        # TODO: I know this is ugly. Suggestions welcome.
        if any(ext in shape_url for ext in archive_types):
            extracted = extract_archive(shape_url, self.workdir)
            for potential_vector in extracted:
                if any(ext in potential_vector for ext in allowed_types):
                    shape_url = potential_vector

        properties = []
        areas = []
        centroid = []

        try:
            with fiona.open(shape_url, 'r', crs=from_epsg(shape_crs)) as src:
                for feature in src:
                    geom = shape(feature['geometry'])

                    properties.append(feature['properties'])
                    centroid.append(geom.centroid.xy)
                    transformed = ops.transform(
                        partial(
                            transform,
                            Proj(src.crs),
                            Proj(from_epsg(projected_crs))),
                        geom)
                    areas.append(transformed.area)

                src.close()
        except Exception as e:
            msg = 'Failed to extract shape from url {}: {}'.format(shape_url, e)
            LOGGER.error(msg)


        response.outputs['properties'] = properties
        response.outputs['area'] = areas
        response['centroid'] = centroid

        return response


# if __name__ == "__main__":
#     from tests.common import TESTDATA
#
#     fields = ['crs={crs}', 'shape=file@xlink:href=file://{file}']
#     datainputs = ';'.join(fields).format(
#         crs=4326,
#         file=TESTDATA['watershed_vector']
#     )
#
#     input = {'crs': 4326, 'file': TESTDATA['watershed_vector']}
#
#     output = {}
#
#     q = ShapeAreaProcess._shapearea_handler(input, output)
#
#     for k in q.keys():
#         print(k, q[k])
