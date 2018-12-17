import logging
from raven.utils import extract_archive
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
                         'Attempt to capture both conatining basin and all upstream basins from point',
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

        lonlat = request.inputs['lonlat_coordinate'][0].data
        try:
            lon, lat = tuple(map(float, re.findall(r'[-+]?[0-9]*\.?[0-9]+', lonlat)))
        except Exception as e:
            msg = 'Failed to parse geo-coordinates {0}: {1}'.format(lonlat, e)
            print(msg)
            logging.error(msg)
            return response

        upstream = request.inputs['collect_upstream'][0].data
        crs = request.inputs['crs'][0].data
        shape_url = request.inputs['shape'][0].file

        archive_types = ['.nc', '.tar', '.zip']
        allowed_types = ['.gml', '.shp', '.geojson', '.json', '.gpkg']

        # TODO: I know this is ugly. Suggestions welcome.
        if any(ext in shape_url for ext in archive_types):
            extracted = extract_archive(shape_url, self.workdir)
            for potential_vector in extracted:
                if any(ext in potential_vector for ext in allowed_types):
                    shape_url = potential_vector

        basin = []
        pfaf = ''
        properties = []
        location = Point(lon, lat)
        try:
            with fiona.open(shape_url, 'r', crs=from_epsg(crs)) as src:
                for feature in src:
                    geometry = shape(feature['geometry'])

                    if geometry.contains(location):
                        basin = feature['properties']['HYBAS_ID']
                        pfaf = feature['properties']['PFAF_ID']
                        prop = {'id': feature['id']}
                        prop.update(feature['properties'])
                        prop.update(feature['geometry'])
                        properties.append(prop)
                        continue

                for feature in src:
                    pfaf_start, pfaf_end = pfaf[0:4], pfaf[4:]
                    basin.append(filter(lambda f: feature['properties']['NEXT_DOWN'] == basin[-1], src))


        except Exception as e:
            msg = 'Failed to extract shape from url {}: {}'.format(shape_url, e)
            LOGGER.error(msg)

        response.outputs['feature'].data = json.dumps(properties)

        return response
