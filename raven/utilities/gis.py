import fiona
from raven.utils import crs_sniffer
from shapely.geometry import shape, Point


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

    Notes
    -----
    This is really slow. Another approach is to use the `fiona.Collection.filter` method.
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


def get_bbox(fn):
    """Return bounding box of first feature in file.

    Parameters
    ----------
    fn : str
      Path to file storing vector features.

    Returns
    -------
    list
      Geographic coordinates of the bounding box (lon0, lat0, lon1, lat1).

    """
    for i, layer_name in enumerate(fiona.listlayers(fn)):
        with fiona.open(fn, 'r', layer=i) as src:
            for feature in src:
                geom = shape(feature['geometry'])
                return geom.bounds

    # for i, layer_name in enumerate(fiona.listlayers(fn)):
    #     with fiona.open(fn, 'r', layer=i) as src:
    #         return src.bounds


def get_dem(bbox, wcs_version='1.0.0'):
    """
    Return a subset of the EarthEnv NorthAmerica DEM.

    Parameters
    ----------
    bbox : sequence
      Geographic coordinates of the bounding box (lon0, lat0, lon1, lat1).
    wcs_version : str
      The OWSlib WCS protocol employed. Default: '1.0.0'

    Returns
    -------
    bytes
      A GeoTIFF array.

    """
    # size: 143999, 87599
    # bbox: -170, 10, -50, 83
    # dlon = dlat = 1200
    from owslib.wcs import WebCoverageService
    from lxml import etree

    crs = 'epsg:4326'
    fmt = 'GeoTIFF'
    (lon0, lat0, lon1, lat1) = bbox

    # We should use 2.0.1. It is supported by GeoServer and the previous versions are deprecated.
    wcs = WebCoverageService('http://boreas.ouranos.ca/geoserver/ows', version=wcs_version)

    # Need to compute the width and height based on the resolution.
    if wcs_version == '1.0.0':
        layer = 'public:EarthEnv_DEM90_NorthAmerica'
        resp = wcs.getCoverage(identifier=layer,
                               bbox=bbox,
                               format=fmt, crs=crs,
                               width=int((lon1 - lon0) * 10), height=int((lat1 - lat0) * 10))

    # This is not working.
    elif wcs_version == '2.0.1':
        layer = "public:EarthEnv_DEM90_NorthAmerica"
        resp = wcs.getCoverage(identifier=[layer, ],
                               format='image/tiff',
                               subsets=[('i', lon0, lon1), ('j', lat0, lat1)])

    else:
        raise NotImplementedError

    data = resp.read()

    try:
        etree.fromstring(data)
        # The response is an XML file describing the server error.
        raise ChildProcessError(data)

    except etree.XMLSyntaxError:
        # The response is the DEM array.
        return data
