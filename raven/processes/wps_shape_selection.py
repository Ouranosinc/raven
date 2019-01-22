import logging
from raven.utils import archive_sniffer
import json
import fiona
import re

from fiona.crs import from_epsg, from_string
from pyproj import Proj, transform
from pywps import LiteralInput, ComplexInput, ComplexOutput
from pywps import Process, FORMATS
from shapely.geometry import shape, Point

LOGGER = logging.getLogger("PYWPS")

LAEA = '+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
WORLDMOLL = '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
ALBERS_NAM = '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs'


class ShapeSelectionProcess(Process):
    """Given lat/lon coordinates and a file containing vector data, return the feature containing the coordinates."""

    def __init__(self):
        inputs = [

            # LiteralInput, coordinates, type string:lat, lon
            # ComplexInput: features, shapefile with multiple features from which one will be selected
            # ComplexOutput: polygonfeature in format that can be easily worked with on the web frontend.

            LiteralInput('collect_upstream',
                         'Attempt to capture both containing basin and all upstream basins from point',
                         data_type='boolean',
                         default='false'),
            LiteralInput('lonlat_coordinate', '(Lon, Lat) tuple for point of interest',
                         data_type='string',
                         default='(-68.724444, 50.646667)'),
            LiteralInput('crs', 'Coordinate Reference System of shape (EPSG) and lat/lon coordinates',
                         data_type='integer',
                         default=4326),
            ComplexInput('shape', 'Vector Shape',
                         abstract='An URL pointing to either an ESRI Shapefile, GML, JSON, GeoJSON.'
                                  ' The ESRI Shapefile must be zipped and contain the .shp, .shx, and .dbf.',
                         min_occurs=1,
                         supported_formats=[FORMATS.GEOJSON, FORMATS.GML, FORMATS.JSON, FORMATS.SHP]),
        ]

        outputs = [
            ComplexOutput('feature', 'Watershed feature geometry',
                          abstract='Geographic representations and descriptions of shape properties.',
                          supported_formats=[FORMATS.GEOJSON]),
            ComplexOutput('upstream_basins', 'IDs for all immediate upstream basins',
                          abstract='List of all first-level tributary basins by their HydroBASINS ID',
                          supported_formats=[FORMATS.JSON])
        ]

        super(ShapeSelectionProcess, self).__init__(
            self._handler,
            identifier="shape-selection",
            title="Shape Selection",
            version="1.0",
            abstract="Return a feature within a polygonal vector file based on lat/lon coordinates.",
            metadata=[],
            inputs=inputs,
            outputs=outputs,
            status_supported=True,
            store_supported=True)

    def _handler(self, request, response):

        # collect_upstream = request['collect_upstream']
        # lonlat = request['lonlat_coordinate']
        # crs = request['crs']
        # shape_url = request['shape']

        collect_upstream = request.inputs['collect_upstream'][0].data
        lonlat = request.inputs['lonlat_coordinate'][0].data
        crs = request.inputs['crs'][0].data
        shape_url = request.inputs['shape'][0].file
        try:
            lon, lat = tuple(map(float, re.findall(r'[-+]?[0-9]*\.?[0-9]+', lonlat)))
        except Exception as e:
            msg = 'Failed to parse geo-coordinates {0}: {1}'.format(lonlat, e)
            LOGGER.error(msg)
            return response

        types = ['.tar', '.zip', '.7z']
        extensions = ['.gml', '.shp', '.geojson', '.json']  # '.gpkg' requires more handling

        shape_url = archive_sniffer(shape_url, working_dir=self.workdir, archive_types=types, extensions=extensions)

        basin = []
        # pfaf = ''
        upstream_basins = []
        properties = []
        location = Point(lon, lat)

        try:
            with fiona.open(shape_url, 'r', crs=from_epsg(crs)) as src:
                for feature in src:
                    geometry = shape(feature['geometry'])

                    if geometry.contains(location):
                        basin = [feature['properties']['HYBAS_ID']]
                        # pfaf = feature['properties']['PFAF_ID']
                        prop = {'id': feature['id']}
                        prop.update(feature['properties'])
                        prop.update(feature['geometry'])
                        properties.append(prop)
                    else:
                        continue

                if collect_upstream:
                    # This can also technically be used to described the drainage network. See HydroBASINS docs.
                    # pfaf_start, pfaf_end = str(pfaf)[0:3], str(pfaf)[3:]

                    up = filter(lambda feature: feature['properties']['NEXT_DOWN'] == basin[0], src)
                    for f in iter(up):
                        upstream_basins.append(f['properties']['HYBAS_ID'])

                src.close()

        except Exception as e:
            msg = 'Failed to extract shape from url {}: {}'.format(shape_url, e)
            LOGGER.error(msg)
            return response

        response.outputs['feature'].data = json.dumps(properties)
        response.outputs['upstream_basins'].data = json.dumps(upstream_basins)
        # response['feature'] = json.dumps(properties)
        # response['upstream_basins'] = upstream_basins

        return response


# if __name__ == '__main__':
#     from tests.common import TESTDATA
#
#     request = {
#         'collect_upstream': True,
#         'lonlat_coordinate': '(-68.724444, 50.646667)',
#         'crs': 4326,
#         'shape': TESTDATA['hydrobasins_12']
#     }
#     response = {}
#
#     s = ShapeSelectionProcess()._handler(request, response)
#     print(response)
