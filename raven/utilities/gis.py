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
    gdf : pd.DataFrame
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

    # Buffer function to fix invalid geometries
    gdf['geometry'] = gdf.buffer(0)

    return gdf.dissolve(by='MAIN_BAS', aggfunc=aggfunc)


def get_bbox(vector, all_features=True):
    """Return bounding box of all features or the first feature in file.

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
        with fiona.open(vector, 'r') as src:
            for feature in src:
                geom = shape(feature['geometry'])
                return geom.bounds

    with fiona.open(vector, 'r') as src:
        return src.bounds


def get_raster_wcs(coordinates, geographic=True, layer=None):
    """Return a subset of a raster image from the local GeoServer via WCS 2.0.1 protocol.

    For geoggraphic rasters, subsetting is based on WGS84 (Long, Lat) boundaries. If not geographic, subsetting based
    on projected coordinate system (Easting, Northing) boundries.

    Parameters
    ----------
    coordinates : sequence
      Geographic coordinates of the bounding box (left, down, right, up)
    geographic : bool
      If True, uses "Long" and "Lat" in WCS call. Otherwise uses "E" and "N".
    layer : str
      Layer name of raster exposed on GeoServer instance. E.g. 'public:CEC_NALCMS_LandUse_2010'

    Returns
    -------
    bytes
      A GeoTIFF array.

    """
    from owslib.wcs import WebCoverageService
    from lxml import etree

    (left, down, right, up) = coordinates

    if geographic:
        x, y = 'Long', 'Lat'
    else:
        x, y = 'E', 'N'

    wcs = WebCoverageService('http://boreas.ouranos.ca/geoserver/ows', version='2.0.1')

    try:
        resp = wcs.getCoverage(identifier=[layer, ],
                               format='image/tiff',
                               subsets=[(x, left, right), (y, down, up)])

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


def get_hydrobasins_location_wfs(coordinates=None, level=12, lakes=True):
    """Return a subset of a raster image from the local GeoServer via WCS 2.0.1 protocol.

    For geoggraphic rasters, subsetting is based on WGS84 (Long, Lat) boundaries. If not geographic, subsetting based
    on projected coordinate system (Easting, Northing) boundries.

    Parameters
    ----------
    coordinates : sequence
      Geographic coordinates of the bounding box (left, down, right, up)
    level : int
      Level of granularity requested for the lakes vector (1:12). Default: 12.
    lakes : bool
      Whether or not the vector should include the delimitation of lakes.

    Returns
    -------
    bytes
      A GeoTIFF array.

    """
    from owslib.wfs import WebFeatureService

    layer = 'public:USGS_HydroBASINS_{}na_lev{}'.format('lake_' if lakes else '', level)

    if coordinates is not None:
        wfs = WebFeatureService('http://boreas.ouranos.ca/geoserver/wfs', version='1.1.0', timeout=30)
        try:
            resp = wfs.getfeature(typename=layer, bbox=coordinates, srsname='urn:x-ogc:def:crs:EPSG:4326')
        except Exception as e:
            raise Exception(e)
    else:
        raise NotImplementedError

    data = resp.read()
    return data


def get_hydrobasins_attributes_wfs(attribute=None, value=None, level=12, lakes=True):
    """Return a subset of a raster image from the local GeoServer via WCS 2.0.1 protocol.

    For geoggraphic rasters, subsetting is based on WGS84 (Long, Lat) boundaries. If not geographic, subsetting based
    on projected coordinate system (Easting, Northing) boundries.

    Parameters
    ----------
    attribute : str
      Attribute/field to be queried.
    value: str
      Value for attribute queried.
    level : int
      Level of granularity requested for the lakes vector (1:12). Default: 12.
    lakes : bool
      Whether or not the vector should include the delimitation of lakes.

    Returns
    -------
    bytes
      A GeoTIFF array.

    """
    from requests import Request
    from owslib.fes import PropertyIsLike
    from lxml import etree

    url = 'http://boreas.ouranos.ca/geoserver/wfs'
    layer = 'public:USGS_HydroBASINS_{}na_lev{}'.format('lake_' if lakes else '', level)

    if attribute is not None and value is not None:

        try:
            attribute = str(attribute)
            value = str(value)

        except ValueError:
            raise Exception('Unable to cast attribute/filter to string')

        try:
            filter_request = PropertyIsLike(propertyname=attribute, literal=value, wildCard='*')
            filterxml = etree.tostring(filter_request.toXML()).decode('utf-8')
            params = dict(service='WFS', version='1.1.0', request='GetFeature', typename=layer, outputFormat='json',
                          filter=filterxml)

            q = Request('GET', url, params=params).prepare().url

        except Exception as e:
            raise Exception(e)

    else:
        raise NotImplementedError

    return q
