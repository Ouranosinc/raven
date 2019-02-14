import json
import logging
import re
import tempfile

import fiona as fio
import geopandas as gpd
from fiona.crs import from_epsg
from pywps import LiteralInput, ComplexInput, ComplexOutput
from pywps import Process, FORMATS
from raven.utils import archive_sniffer, crs_sniffer
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
            ComplexOutput('geojson', 'Watershed feature geometry',
                          abstract='Geographic representations of shape properties.',
                          supported_formats=[FORMATS.GEOJSON]),
            ComplexOutput('upstream_basins', 'HydroBASINS IDs for all immediate upstream basins',
                          abstract='Exhaustive list of all tributary sub-basins according to their HydroBASINS IDs',
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

        def get_upstream_basins(basin_dataframe, basin_id):
            up_str_mask = basin_dataframe['NEXT_DOWN'] == basin_id
            iter_upstream = basin_dataframe[up_str_mask]['HYBAS_ID']
            return iter_upstream

        collect_upstream = request['collect_upstream']
        lonlat = request['lonlat_coordinate']
        crs = from_epsg(request['crs'])
        shape_url = request['shape']

        # collect_upstream = request.inputs['collect_upstream'][0].data
        # lonlat = request.inputs['lonlat_coordinate'][0].data
        # crs = from_epsg(request.inputs['crs'][0].data)
        # shape_url = request.inputs['shape'][0].file

        try:
            lon, lat = tuple(map(float, re.findall(r'[-+]?[0-9]*\.?[0-9]+', lonlat)))
        except Exception as e:
            msg = 'Failed to parse geo-coordinates {}: {}'.format(lonlat, e)
            LOGGER.error(msg)
            return response

        extensions = ['.gml', '.shp', '.geojson', '.json']  # '.gpkg' requires more handling
        shape_url = archive_sniffer(shape_url, working_dir=self.workdir, extensions=extensions)

        basin = []
        geojson = tempfile.NamedTemporaryFile(suffix='.json')
        upstream_basins = []
        location = Point(lon, lat)
        crs = crs_sniffer(shape_url, crs)

        with fio.Env():  # hacky workaround for pip-installed fiona; Can be removed in conda systems
            try:
                with fio.open(shape_url, 'r', crs=crs) as src:
                    for feat in iter(src):
                        geom = shape(feat['geometry'])

                        if geom.contains(location):
                            # pfaf = feature['properties']['PFAF_ID']
                            if collect_upstream:
                                basin = feat['properties']['HYBAS_ID']
                            else:
                                with open(geojson, 'w') as f:
                                    json.dump(feat, f)
                            continue
                    src.close()

                if collect_upstream:
                    with fio.Collection(shape_url, 'r', crs=crs) as src:
                        gdf = gpd.GeoDataFrame.from_features(src, crs=crs)
                        row_id = gdf['HYBAS_ID'] == basin
                        main_basin_id = gdf[row_id]['MAIN_BAS']
                        main_basin_mask = gdf['MAIN_BAS'] == main_basin_id.values[0]
                        main_basin_gdf = gdf[main_basin_mask]

                        all_basins = list(main_basin_gdf[row_id]['HYBAS_ID'].values)
                        for i in all_basins:
                            tmp = get_upstream_basins(main_basin_gdf, i)
                            if len(tmp):
                                all_basins.extend(tmp.values)

                        df_sub = main_basin_gdf[main_basin_gdf['HYBAS_ID'] == all_basins[0]]
                        for a in all_basins[1:]:
                            df_sub = df_sub.append(main_basin_gdf[main_basin_gdf['HYBAS_ID'] == a])
                        dissolved = df_sub.dissolve(by='MAIN_BAS').to_json()
                        upstream_basins = all_basins

                        LOGGER.warning('Writing to {}'.format(geojson.name))
                        with open(geojson.name, 'w') as f:
                            json.dump(dissolved, f) # TODO: Figure out why this doesn't work

                        # fn = os.path.join(self.workdir, 'upstream.json')
                        # with open(fn, 'w') as f:
                        #     json.dump(stats, f)
                        # response.outputs['properties'].file = fn

                        # filename = '{}.shp'.format(basin)

                        # dissolved.to_file('output.json', driver='GeoJSON')

                    # with fio.open(shape_url, 'r', crs=from_epsg(crs)) as src:
                    #     # 'PFAF_ID' can also technically be used to described the drainage network. See HydroBASINS docs.
                    #     # pfaf_start, pfaf_end = str(pfaf)[0:3], str(pfaf)[3:]
                    #     print(len(src))
                    #     upstream = filter(lambda f: feature['properties']['NEXT_DOWN'] == basin[0], src)

            except Exception as e:
                msg = 'Failed to extract shape from url {}: {}'.format(shape_url, e)
                LOGGER.error(msg)
                return response

        # response.outputs['geojson'].file = geojson
        # response.outputs['upstream_basins'].data = upstream_basins

        response['geojson'] = geojson
        response['upstream_basins'] = upstream_basins

        return response


def testing():
    inputs = dict(collect_upstream=True,
                  lonlat_coordinate="(-68.724444, 50.646667)",
                  crs=4326,
                  shape='/home/tjs/git/raven/tests/testdata/usgs_hydrobasins/hybas_lake_na_lev12_v1c.zip')
    outputs = {}
    ShapeSelectionProcess()._handler(request=inputs, response=outputs)

    print(outputs['geojson'], outputs['upstream_basins'])


if __name__ == '__main__':
    testing()


