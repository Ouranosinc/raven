import json
import logging
import re

import fiona as fio
from fiona.crs import from_epsg
from pywps import LiteralInput, ComplexInput, ComplexOutput
from pywps import Process, FORMATS
from raven.utils import archive_sniffer
from shapely.geometry import shape, Point

LOGGER = logging.getLogger("PYWPS")


class ShapeSelectionProcess(Process):
    """Given lat/lon coordinates and a file containing vector data, return the feature containing the coordinates."""

    def __init__(self):
        inputs = [
            LiteralInput('lonlat_coordinate', '(Lon, Lat) tuple for point of interest',
                         data_type='string',
                         default='(-68.724444, 50.646667)'),
            LiteralInput('collect_upstream',
                         'Attempt to capture both containing basin and all upstream basins from point',
                         data_type='boolean',
                         default='false'),
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
        upstream_basins = []
        properties = []
        location = Point(lon, lat)

        try:
            with fio.open(shape_url, 'r', crs=from_epsg(crs)) as src:
                for feature in src:
                    geom = shape(feature['geometry'])

                    if geom.contains(location):
                        # pfaf = feature['properties']['PFAF_ID']

                        basin = [feature['properties']['HYBAS_ID']]
                        prop = {'id': feature['id']}
                        prop.update(feature['properties'])
                        prop.update(feature['geometry'])
                        properties.append(prop)

                    else:
                        continue

                if collect_upstream:
                    # 'PFAF_ID' can also technically be used to described the drainage network. See HydroBASINS docs.
                    # pfaf_start, pfaf_end = str(pfaf)[0:3], str(pfaf)[3:]

                    upstream = filter(lambda f: feature['properties']['NEXT_DOWN'] == basin[0], src)
                    for up in iter(upstream):
                        upstream_basins.append(up['properties']['HYBAS_ID'])

                src.close()

        except Exception as e:
            msg = 'Failed to extract shape from url {}: {}'.format(shape_url, e)
            LOGGER.error(msg)
            return response

        # response['feature'] = json.dumps(properties)
        # response['upstream_basins'] = upstream_basins
        response.outputs['feature'].data = json.dumps(properties)
        response.outputs['upstream_basins'].data = json.dumps(upstream_basins)

        return response
