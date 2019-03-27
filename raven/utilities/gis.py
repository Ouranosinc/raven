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
    dict
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


def hydrobasins_upstream_ids(fid, df):
    """Return a list of hydrobasins features located upstream.

    Parameters
    ----------
    fid : feature id
      HYBAS_ID of the downstream feature.
    df : pd.DataFrame
      Watershed attributes.

    Returns
    -------
    pd.Series
      Basins ids including `fid` and its upstream contributors.
    """

    def upstream_ids(bdf, bid):
        return bdf[bdf['NEXT_DOWN'] == bid]['HYBAS_ID']

    # Locate the downstream feature
    ds = df.set_index("HYBAS_ID").loc[fid]

    # Do a first selection on the main basin ID of the downstream feature.
    sub = df[df['MAIN_BAS'] == ds['MAIN_BAS']]

    # Find upstream basins
    up = [fid, ]
    for b in up:
        tmp = upstream_ids(sub, b)
        if len(tmp):
            up.extend(tmp)

    return sub[sub['HYBAS_ID'].isin(up)]


def hydrobasins_aggregate(gdf):
    """Aggregate multiple hydrobasin watersheds into a single geometry.

    Parameters
    ----------
    ids : sequence
      Basins ids, namely the HYBAS_ID attribute.
    df : pd.DataFrame
      Watershed attributes indexed by HYBAS_ID

    Returns
    -------

    """
    # TODO: Review. Not sure it all makes sense.
    def aggfunc(x):
        if x.name in ['COAST', 'DIST_MAIN', 'DIST_SINK']:
            return x.min()
        elif x.name in ['SUB_AREA', 'LAKE']:
            return x.sum()
        else:
            return x[0]

    return gdf.dissolve(by='MAIN_BAS', aggfunc=aggfunc)
