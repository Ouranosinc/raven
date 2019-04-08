# -*- coding: utf-8 -*-
import json
import logging
from pathlib import Path

from pywps import LiteralInput, ComplexOutput
from pywps import Process, FORMATS
from raven.utils import archive_sniffer, single_file_check, parse_lonlat
from shapely.geometry import Point
from raven.utilities import gis
from raven.utils import crs_sniffer
import fiona
import pandas as pd
import geopandas as gpd

from tests.common import TESTDATA

LOGGER = logging.getLogger("PYWPS")

DATA = Path(__file__).parent / 'hydrobasins_tables'
HYBAS = {
    'lake_lev07': DATA / 'hybas_lake_na_lev07.csv',
    'lake_lev08': DATA / 'hybas_lake_na_lev08.csv',
    'lake_lev09': DATA / 'hybas_lake_na_lev09.csv',
    'lake_lev10': DATA / 'hybas_lake_na_lev10.csv',
    'lake_lev11': DATA / 'hybas_lake_na_lev11.csv',
    'lake_lev12': DATA / 'hybas_lake_na_lev12.csv',
    'lev07': DATA / 'hybas_na_lev07.csv',
    'lev08': DATA / 'hybas_na_lev08.csv',
    'lev09': DATA / 'hybas_na_lev09.csv',
    'lev10': DATA / 'hybas_na_lev10.csv',
    'lev11': DATA / 'hybas_na_lev11.csv',
    'lev12': DATA / 'hybas_na_lev12.csv',
}


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

        shape_description = 'hydrobasins_{}na_lev{}'.format('lake_' if lakes else '', level)
        table = DATA / 'hybas_{}na_lev{:02}.csv'.format('lake_' if lakes else '', level)
        shape_url = TESTDATA[shape_description]

        # Why not just map(float, lonlat.split(',')) ?
        lon, lat = parse_lonlat(lonlat)
        bbox = (lon, lat, lon, lat)

        extensions = ['.gml', '.shp', '.gpkg', '.geojson', '.json']
        shp = single_file_check(archive_sniffer(shape_url, working_dir=self.workdir, extensions=extensions))

        shape_crs = crs_sniffer(shp)
        with fiona.Collection(shp, 'r', crs=shape_crs) as src:

            # Find feature containing location
            feat = next(src.filter(bbox=bbox))
            hybas_id = feat['properties']['HYBAS_ID']

            # check conditions
            if collect_upstream:
                if lakes is False or level != 12:
                    raise ValueError("Set lakes to True and level to 12.")

                # Read table of attributes relevant to upstream watershed identification
                df = pd.read_csv(table)

                # Get upstream features, including feat
                up = gis.hydrobasins_upstream_ids(hybas_id, df)

                # Read only upstream geometries
                gdf = gpd.GeoDataFrame.from_features((src.get(i) for i in up.index))

                # Aggregate upstream features into a single geometry.
                agg = gis.hydrobasins_aggregate(gdf)

                response.outputs['feature'].data = agg.to_json()
                response.outputs['upstream_ids'].data = json.dumps(up['HYBAS_ID'].tolist())

            else:
                response.outputs['feature'].data = json.dumps(feat)
                response.outputs['upstream_ids'].data = json.dumps([hybas_id, ])

        return response
