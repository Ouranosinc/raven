# -*- coding: utf-8 -*-
import json
import logging
import tempfile

import fiona as fio
import geopandas as gpd

from pywps import LiteralInput, ComplexOutput
from pywps import Process, FORMATS
from raven.utils import archive_sniffer, crs_sniffer, single_file_check, parse_lonlat
from shapely.geometry import shape, Point

from tests.common import TESTDATA

LOGGER = logging.getLogger("PYWPS")


class ShapeSelectionProcess(Process):
    """Given lat/lon coordinates and a file containing vector data, return the feature containing the coordinates."""

    def __init__(self):
        inputs = [
            LiteralInput('lonlat_coordinate', '(Lon, Lat) tuple for point of interest',
                         data_type='string',
                         default='(-68.724444, 50.646667)',
                         min_occurs=1, max_occurs=1),
            LiteralInput('level', 'Resolution level of HydroBASINS Shapes',
                         data_type='integer',
                         default=12,
                         allowed_values=[7, 8, 9, 10, 11, 12],
                         min_occurs=1, max_occurs=1),
            LiteralInput('lakes', 'Use the HydroBASINS version that includes lake outlines',
                         data_type='boolean',
                         default='true',
                         min_occurs=1, max_occurs=1),
            LiteralInput('collect_upstream',
                         'Attempt to capture both the containing basin and all tributary basins from point',
                         data_type='boolean',
                         default='false',
                         min_occurs=1, max_occurs=1),
        ]

        outputs = [
            ComplexOutput('geojson', 'Watershed feature geometry',
                          abstract='Geographic representations of shape properties.',
                          supported_formats=[FORMATS.JSON]),
            ComplexOutput('upstream_basins', 'HydroBASINS IDs for all immediate upstream basins',
                          abstract='Exhaustive list of all tributary sub-basins according to their HydroBASINS IDs',
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

        def get_upstream_hydrobasins(basin_dataframe, basin_id):
            up_str_mask = basin_dataframe['NEXT_DOWN'] == basin_id
            iter_upstream = basin_dataframe[up_str_mask]['HYBAS_ID']
            return iter_upstream

        level = request.inputs['level'][0].data
        lakes = request.inputs['lakes'][0].data
        collect_upstream = request.inputs['collect_upstream'][0].data
        lonlat = request.inputs['lonlat_coordinate'][0].data
        lon, lat = parse_lonlat(lonlat)

        if lakes:
            shape_description = 'hydrobasins_lake_na_lev{}'.format(level)
        else:
            shape_description = 'hydrobasins_na_lev{}'.format(level)

        if shape_description != 'hydrobasins_lake_na_lev12':
            msg = 'Options for HydroBASINS levels and lakes will be in prod'
            LOGGER.error(NotImplementedError(msg))
            shape_description = 'hydrobasins_lake_na_lev12'

        shape_url = TESTDATA[shape_description]

        extensions = ['.gml', '.shp', '.geojson', '.json']  # '.gpkg' requires more handling
        vector_file = single_file_check(archive_sniffer(shape_url, working_dir=self.workdir, extensions=extensions))

        basin = []
        upstream_basins = []
        location = Point(lon, lat)
        found = False

        shape_crs = crs_sniffer(vector_file)

        with fio.Env():  # Workaround for pip-installed fiona; Can be removed in conda-installed systems
            geojson = tempfile.NamedTemporaryFile(suffix='.json', delete=False)
            try:
                with fio.open(vector_file, 'r', crs=shape_crs) as src:
                    for feat in iter(src):
                        geom = shape(feat['geometry'])

                        if geom.contains(location):
                            found = True
                            # pfaf = feature['properties']['PFAF_ID']
                            if collect_upstream:
                                basin = feat['properties']['HYBAS_ID']
                            else:
                                LOGGER.info('Writing feature to {}'.format(geojson.name))
                                with open(geojson.name, 'w') as f:
                                    json.dump(feat, f)
                            break
                    src.close()
                if not found:
                    msg = 'No basin found at lon:{}, lat:{}'.format(lon, lat)
                    LOGGER.exception(msg)
                    raise ValueError(msg)

                if collect_upstream:
                    LOGGER.info('Collecting upstream from basin {}'.format(basin))
                    with fio.Collection(vector_file, 'r', crs=shape_crs) as src:
                        gdf = gpd.GeoDataFrame.from_features(src, crs=shape_crs)
                        row_id = gdf['HYBAS_ID'] == basin
                        main_basin_id = gdf[row_id]['MAIN_BAS']
                        main_basin_mask = gdf['MAIN_BAS'] == main_basin_id.values[0]
                        main_basin_gdf = gdf[main_basin_mask]

                        all_basins = list(main_basin_gdf[row_id]['HYBAS_ID'].values)
                        for i in all_basins:
                            tmp = get_upstream_hydrobasins(main_basin_gdf, i)
                            if len(tmp):
                                all_basins.extend(tmp.values)

                        upstream_basins = [x.item() for x in all_basins]  # Convert from numpy Int64

                        df_sub = main_basin_gdf[main_basin_gdf['HYBAS_ID'] == all_basins[0]]
                        for a in all_basins[1:]:
                            df_sub = df_sub.append(main_basin_gdf[main_basin_gdf['HYBAS_ID'] == a])
                        dissolved = df_sub.dissolve(by='MAIN_BAS').to_json()

                        LOGGER.info('Writing union of features to {}'.format(geojson.name))
                        with open(geojson.name, 'w') as f:
                            f.write(dissolved)

                    # with fio.open(shape_url, 'r', crs=from_epsg(crs)) as src:
                    #     # 'PFAF_ID' can also technically be used to described the drainage network
                    #       See HydroBASINS docs.
                    #     # pfaf_start, pfaf_end = str(pfaf)[0:3], str(pfaf)[3:]

            except Exception as e:
                msg = '{}: Failed to perform analysis using {} and location {}'.format(e, shape_url, lon, lat)
                LOGGER.error(msg)
                raise Exception(msg)

        response.outputs['geojson'].data = geojson.name
        response.outputs['upstream_basins'].data = json.dumps(upstream_basins)

        return response
