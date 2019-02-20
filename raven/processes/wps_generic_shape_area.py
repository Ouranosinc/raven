import fiona
import json
import logging
from fiona.crs import from_epsg
from rasterio.crs import CRS
from pywps import LiteralInput, ComplexInput, ComplexOutput
from pywps import Process, FORMATS
from raven.utils import archive_sniffer, geom_transform, geom_prop, equal_area_geom_prop
from shapely.geometry import shape

LOGGER = logging.getLogger("PYWPS")


class ShapeAreaProcess(Process):
    """Given a file containing vector data, provide general information and spatial characteristics"""

    def __init__(self):
        inputs = [
            ComplexInput('shape', 'Vector Shape',
                         abstract='An ESRI Shapefile, GML, JSON or GeoJSON file. The ESRI Shapefile must be zipped and'
                                  ' contain the .shp, .shx, and .dbf.',
                         supported_formats=[FORMATS.GML, FORMATS.GEOJSON, FORMATS.SHP, FORMATS.JSON],
                         min_occurs=1, max_occurs=1),
            LiteralInput('crs', 'Coordinate Reference System of shape (EPSG code; Default: 4326)',
                         data_type='integer',
                         default=4326,
                         min_occurs=1, max_occurs=1),
            LiteralInput('projected_crs',
                         'Coordinate Reference System for area calculation (EPSG code; Default:32198)',
                         data_type='integer',
                         default=32198,
                         min_occurs=1, max_occurs=1)
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

        extensions = ['.gml', '.shp', '.geojson', '.json']  # '.gpkg' requires more handling

        shape_url = archive_sniffer(shape_url, working_dir=self.workdir, extensions=extensions)

        try:
            projection = CRS.from_user_input(projected_crs)
            if not projection.is_projected:
                msg = 'Destination CRS {} is not projected.' \
                      'Terrain analysis values may be erroneous.'.format(projection.to_epsg())
                LOGGER.warning(msg)
        except Exception as e:
            msg = '{}: Failed to parse CRS definition. Exiting.'.format(e)
            LOGGER.error(msg)
            return response

        properties = []

        try:
            with fiona.open(shape_url, 'r', crs=from_epsg(shape_crs)) as src:
                for feature in src:
                    geom = shape(feature['geometry'])

                    transformed = geom_transform(geom, source_crs=shape_crs, target_crs=projected_crs)
                    prop = {'id': feature['id']}
                    prop.update(feature['properties'])
                    prop.update(geom_prop(geom))
                    prop.update(equal_area_geom_prop(transformed))
                    properties.append(prop)
                    break

        except Exception as e:
            msg = 'Failed to extract shape from url {}: {}'.format(shape_url, e)
            LOGGER.error(msg)

        response.outputs['properties'].data = json.dumps(properties)

        return response
