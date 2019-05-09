import fiona
import json
import logging

from pywps import LiteralInput, ComplexInput, ComplexOutput
from pywps import Process, FORMATS
from raven.utils import archive_sniffer, crs_sniffer, single_file_check
from raven.utils import geom_transform, geom_prop
from rasterio.crs import CRS
from shapely.geometry import shape

LOGGER = logging.getLogger("PYWPS")


class ShapePropertiesProcess(Process):
    """Given a file containing vector data, provide general information and spatial characteristics."""

    def __init__(self):
        inputs = [
            ComplexInput('shape', 'Vector Shape',
                         abstract='An ESRI Shapefile, GML, GeoPackage, JSON or GeoJSON file.'
                                  ' The ESRI Shapefile must be zipped and contain the .shp, .shx, and .dbf.',
                         supported_formats=[FORMATS.GML, FORMATS.GEOJSON, FORMATS.SHP, FORMATS.JSON],
                         min_occurs=1, max_occurs=1),
            LiteralInput('projected_crs',
                         'Coordinate Reference System for area calculation (Default: EPSG:6622,'
                         ' NAD83(CSRS) / Quebec Lambert)',
                         data_type='integer',
                         default=6622,
                         min_occurs=1, max_occurs=1)]
        outputs = [
            ComplexOutput('properties', 'Feature schemas',
                          abstract='Geographic representations and descriptions of shape properties: '
                                   'centroid coordinates, area, perimeter and gravelius index.',
                          supported_formats=[FORMATS.JSON],
                          ),
        ]

        super(ShapePropertiesProcess, self).__init__(
            self._handler,
            identifier="shape-properties",
            title="Shape Properties",
            version="1.0",
            abstract="Return shape area in square metres based on line boundaries of a polygonal vector file.",
            metadata=[],
            inputs=inputs,
            outputs=outputs,
            status_supported=True,
            store_supported=True)

    def _handler(self, request, response):

        shape_url = request.inputs['shape'][0].file
        projected_crs = request.inputs['projected_crs'][0].data

        extensions = ['.gml', '.shp', '.gpkg', '.geojson', '.json']
        vector_file = single_file_check(archive_sniffer(shape_url, working_dir=self.workdir, extensions=extensions))
        shape_crs = crs_sniffer(vector_file)

        try:
            projection = CRS.from_epsg(projected_crs)
            if projection.is_geographic:
                msg = 'Desired CRS {} is geographic. ' \
                      'Areal analysis values will be in decimal-degree units.'.format(projection.to_epsg())
                LOGGER.warning(msg)
        except Exception as e:
            msg = '{}: Failed to parse CRS definition. Exiting.'.format(e)
            LOGGER.error(msg)
            raise Exception(msg)

        properties = []
        try:
            for i, layer_name in enumerate(fiona.listlayers(vector_file)):
                with fiona.open(vector_file, 'r', crs=shape_crs, layer=i) as src:
                    for feature in src:
                        geom = shape(feature['geometry'])

                        transformed = geom_transform(geom, source_crs=shape_crs, target_crs=projection)
                        prop = {'id': feature['id']}
                        prop.update(feature['properties'])
                        prop.update(geom_prop(transformed))

                        # Recompute the centroid location using the original projection
                        prop['centroid'] = geom_prop(geom)['centroid']

                        properties.append(prop)

        except Exception as e:
            msg = '{}: Failed to extract features from shape {}'.format(e, vector_file)
            LOGGER.error(msg)
            raise Exception(msg)

        response.outputs['properties'].data = json.dumps(properties)

        return response
