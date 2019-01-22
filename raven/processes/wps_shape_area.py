import fiona
import json
import logging
from fiona.crs import from_epsg
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
                         abstract='An URL pointing to either an ESRI Shapefile, GML, JSON, GeoJSON.'
                                  ' The ESRI Shapefile must be zipped and contain the .shp, .shx, and .dbf.',
                         min_occurs=1,
                         supported_formats=[FORMATS.GML, FORMATS.GEOJSON, FORMATS.SHP, FORMATS.JSON]),
            LiteralInput('crs', 'Coordinate Reference System of shape (EPSG code; Default: 4326)',
                         data_type='integer',
                         default=4326),
            LiteralInput('projected_crs',
                         'Coordinate Reference System for area calculation (EPSG code; Default:32198)',
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

        # shape_url = request['shape']
        # shape_crs = request['crs']
        # projected_crs = request['projected_crs']

        shape_url = request.inputs['shape'][0].file
        shape_crs = request.inputs['crs'][0].data
        projected_crs = request.inputs['projected_crs'][0].data

        types = ['.tar', '.zip', '.7z']
        extensions = ['.gml', '.shp', '.geojson', '.json']  # '.gpkg' requires more handling

        shape_url = archive_sniffer(shape_url, working_dir=self.workdir, archive_types=types, extensions=extensions)

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

        # response['properties'] = properties
        response.outputs['properties'].data = json.dumps(properties)

        return response


# if __name__ == '__main__':
#     from tests.common import TESTDATA
#
#     request = {
#         'shape': TESTDATA['hydrobasins_12'],
#         'crs': 4326,
#         'projected_crs': 32198
#     }
#     response = {}
#
#     s = ShapeAreaProcess()._handler(request, response)
#     print(response['properties'])
