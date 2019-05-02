import fiona
import collections
from raven.utils import crs_sniffer, single_file_check
from shapely.geometry import shape, Point


def feature_contains(point, shp):
    """Return the first feature containing a location.

    Parameters
    ----------
    point : shapely.Point or tuple
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

    if isinstance(point, collections.abc.Sequence) and not isinstance(point, str):
        for coord in point:
            if isinstance(coord, (int, float)):
                pass
    elif isinstance(point, Point):
        pass
    else:
        raise ValueError("point should be shapely.Point or tuple of coordinates, got : {}".format(point))

    shape_crs = crs_sniffer(single_file_check(shp))

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


def get_bbox(vector, all_features=True):
    """Return bounding box of first feature in file.

    Parameters
    ----------
    vector : str
      Path to file storing vector features.
    all_features : bool
      Return the bounding box for all features. Default: True.

    Returns
    -------
    list
      Geographic coordinates of the bounding box (lon0, lat0, lon1, lat1).

    """

    if not all_features:
        for i, layer_name in enumerate(fiona.listlayers(vector)):
            with fiona.open(vector, 'r', layer=i) as src:
                for feature in src:
                    geom = shape(feature['geometry'])
                    return geom.bounds

    for i, layer_name in enumerate(fiona.listlayers(vector)):
        with fiona.open(vector, 'r', layer=i) as src:
            return src.bounds


def get_dem_wcs(bbox):
    """Return a subset of the EarthEnv NorthAmerica DEM.

    Subsetting based on WGS84 (Long, Lat) boundaries

    Parameters
    ----------
    bbox : sequence
      Geographic coordinates of the bounding box (lon0, lat0, lon1, lat1).

    Returns
    -------
    bytes
      A GeoTIFF array.

    """
    from owslib.wcs import WebCoverageService
    from lxml import etree

    (lon0, lat0, lon1, lat1) = bbox

    wcs = WebCoverageService('http://boreas.ouranos.ca/geoserver/ows', version='2.0.1')

    try:
        layer = "public:EarthEnv_DEM90_NorthAmerica"
        resp = wcs.getCoverage(identifier=[layer, ],
                               format='image/tiff',
                               subsets=[('Long', lon0, lon1), ('Lat', lat0, lat1)])

    except Exception as e:
        raise Exception(e)

    data = resp.read()

    try:
        etree.fromstring(data)
        # The response is an XML file describing the server error.
        raise ChildProcessError(data)

    except etree.XMLSyntaxError:
        # The response is the DEM array.
        return data


def get_nalcms_wcs(bbox, year=2010):
    """Return a subset of the CEC North American Land Change Monitoring System image.

    Subsetting based on(Easting, Northing) using a LAEA projected coordinate system.

    Parameters
    ----------
    bbox : sequence
      Geographic coordinates of the bounding box (E0, N0, E1, N1).
    year : int
      Year of interest. Default: 2010.

    Returns
    -------
    bytes
      A GeoTIFF array.

    """
    from owslib.wcs import WebCoverageService
    from lxml import etree

    (E0, N0, E1, N1) = bbox

    wcs = WebCoverageService('http://boreas.ouranos.ca/geoserver/ows', version='2.0.1')

    if year not in [2005, 2010]:
        raise NotImplementedError

    try:
        layer = "public:CEC_NALCMS_LandUse_{}".format(year)
        resp = wcs.getCoverage(identifier=[layer, ],
                               format='image/tiff',
                               subsets=[('E', E0, E1), ('N', N0, N1)])

    except Exception as e:
        raise Exception(e)

    data = resp.read()

    try:
        etree.fromstring(data)
        # The response is an XML file describing the server error.
        raise ChildProcessError(data)

    except etree.XMLSyntaxError:
        # The response is the DEM array.
        return data
