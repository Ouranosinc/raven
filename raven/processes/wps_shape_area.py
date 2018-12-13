import logging
from functools import partial
from raven.utils import extract_archive, geom_prop, equal_area_geom_prop
import json
import fiona
import shapely.ops as ops
from fiona.crs import from_epsg, from_string
from pyproj import Proj, transform
from pywps import LiteralInput, ComplexInput, ComplexOutput
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
            ComplexOutput('properties', 'Feature schemas',
                          abstract='Geographic representations and descriptions of shape properties.',
                          supported_formats=[FORMATS.JSON],
                          ),
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

        try:
            with fiona.open(shape_url, 'r', crs=from_epsg(shape_crs)) as src:
                for feature in src:

                    geom = shape(feature['geometry'])
                    transformed = ops.transform(
                        partial(
                            transform,
                            Proj(src.crs),
                            Proj(from_epsg(projected_crs))),
                        geom)

                    prop = {'id': feature['id']}
                    prop.update(feature['properties'])
                    prop.update(geom_prop(geom))
                    prop.update(equal_area_geom_prop(transformed))
                    properties.append(prop)

        except Exception as e:
            msg = 'Failed to extract shape from url {}: {}'.format(shape_url, e)
            LOGGER.error(msg)

        response.outputs['properties'].data = json.dumps(properties)

        return response
