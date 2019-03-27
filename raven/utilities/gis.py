import fiona
from raven.utils import crs_sniffer
from shapely.geometry import shape, Point
import geopandas as gpd


def feature_contains(point, shp):
    """Return the first feature containing a location.

    Parameters
    ----------
    point : shapely.Point
      Location coordinates.
    shp : str
      Path to the file storing the geometries.

    Returns
    -------
    str
      The feature found.
    """
    if not isinstance(point, Point):
        raise ValueError("point should be shapely.Point instance, got : {}".format(point))

    shape_crs = crs_sniffer(shp)
    with fiona.Env():
        for i, layer_name in enumerate(fiona.listlayers(shp)):
            with fiona.open(shp, 'r', crs=shape_crs, layer=i) as src:
                for feat in iter(src):
                    geom = shape(feat['geometry'])

                    if geom.contains(point):
                        return feat

    raise LookupError("Could not find feature containing point {} in {}.".format(point, shp))


def hydrobasins_upstream_features(fid, df):
    """Return a list of hydrobasins features located upstream.

    Parameters
    ----------
    fid : feature id
      HYBAS_ID of the downstream feature.
    df : pd.DataFrame
      Watershed attributes.

    Returns
    -------
    list
      Basins ids including `fid` and its upstream contributors.
    """

    df.set_index('HYBAS_ID', inplace=True)

    def upstream_id(bdf, bid):
        return bdf[bdf['NEXT_DOWN'] == bid].index.values.tolist()

    # Locate the downstream feature
    ds = df.loc[fid]

    # Do a first selection on the main basin ID of the downstream feature.
    sub = df[df['MAIN_BAS'] == ds['MAIN_BAS']]

    # Find upstream basins
    up = [fid, ]
    for b in up:
        tmp = upstream_id(sub, b)
        if len(tmp):
            up.extend(tmp)

    return up


def hydrobasins_aggregate(ids, shp):
    """Aggregate multiple hydrobasin watersheds into a single geometry.

    Parameters
    ----------
    ids : sequence
      Basins ids, namely the HYBAS_ID attribute.
    shp : path
      Path to shapefile storing the hydrobasins geometries.

    Returns
    -------

    """

    shape_crs = crs_sniffer(shp)
    with fiona.Collection(shp, 'r', crs=shape_crs) as src:
        gdf = gpd.GeoDataFrame.from_features(src, crs=shape_crs).set_index('HYBAS_ID')
        up = gdf.loc[ids]
        return up.dissolve(by='MAIN_BAS', aggfunc='sum')
