# -*- coding: utf-8 -*-
import json
import logging

from pywps import LiteralInput, ComplexOutput
from pywps import Process, FORMATS
from raven.utils import archive_sniffer, single_file_check, parse_lonlat
from shapely.geometry import Point
from raven.utilities import gis

from tests.common import TESTDATA

LOGGER = logging.getLogger("PYWPS")


class ShapeSelectionProcess(Process):
    """Given lat/lon coordinates and a file containing vector data, return the feature containing the coordinates."""

    def __init__(self):
        inputs = [
            LiteralInput('location', 'Location coordinates (lon, lat)',
                         abstract="Location coordinates (longitude, latitude) for point of interest.",
                         data_type='string',
                         default='-68.724444, 50.646667',
                         min_occurs=1,
                         max_occurs=1),
            LiteralInput('level', 'Resolution level of HydroBASINS Shapes',
                         data_type='integer',
                         default=12,
                         allowed_values=[7, 8, 9, 10, 11, 12],
                         min_occurs=1,
                         max_occurs=1),
            LiteralInput('lakes', 'Use the HydroBASINS version that includes lake outlines',
                         data_type='boolean',
                         default='true',
                         min_occurs=1,
                         max_occurs=1),
            LiteralInput('aggregate_upstream',
                         'Attempt to capture both the containing basin and all tributary basins from point',
                         data_type='boolean',
                         default='false',
                         min_occurs=1,
                         max_occurs=1),
        ]

        outputs = [
            ComplexOutput('feature', 'Watershed feature geometry',
                          abstract='Geographic representation of shape properties.',
                          supported_formats=[FORMATS.GEOJSON]),
            ComplexOutput('upstream_ids', 'HydroBASINS IDs for all immediate upstream basins',
                          abstract='List of all tributary sub-basins according to their HydroBASINS IDs, '
                                   'including the downstream basin.',
                          supported_formats=[FORMATS.JSON])
        ]

        super(ShapeSelectionProcess, self).__init__(
            self._handler,
            identifier="shape-selection",
            title="Shape Selection",
            version="1.0",
            abstract="Return a feature within a polygon vector file based on lon/lat coordinates.",
            metadata=[],
            inputs=inputs,
            outputs=outputs,
            status_supported=True,
            store_supported=True)

    def _handler(self, request, response):

        level = request.inputs['level'][0].data
        lakes = request.inputs['lakes'][0].data
        collect_upstream = request.inputs['aggregate_upstream'][0].data
        lonlat = request.inputs['location'][0].data

        if lakes:
            shape_description = 'hydrobasins_lake_na_lev{}'.format(level)
        else:
            shape_description = 'hydrobasins_na_lev{}'.format(level)

        if shape_description != 'hydrobasins_lake_na_lev12':
            msg = 'Options for HydroBASINS levels and lakes will be in prod'
            LOGGER.error(NotImplementedError(msg))
            shape_description = 'hydrobasins_lake_na_lev12'

        shape_url = TESTDATA[shape_description]

        # Why not just map(float, lonlat.split(',')) ?
        lon, lat = parse_lonlat(lonlat)
        location = Point(lon, lat)

        extensions = ['.gml', '.shp', '.gpkg', '.geojson', '.json']
        vector_file = single_file_check(archive_sniffer(shape_url, working_dir=self.workdir, extensions=extensions))

        # Find feature containing location
        feat = gis.feature_contains(location, vector_file)
        hybas_id = feat['properties']['HYBAS_ID']

        if collect_upstream:
            # Get upstream features, including feat
            up = gis.hydrobasins_upstream_features(hybas_id, vector_file)

            # Aggregate upstream features into a single geometry.
            response.outputs['feature'].data = up.dissolve(by='MAIN_BAS').to_json()
            response.outputs['upstream_ids'].data = json.dumps(up.index.values.tolist())
        else:
            response.outputs['feature'].data = json.dumps(feat)
            response.outputs['upstream_ids'].data = json.dumps([hybas_id, ])

        return response
