import fiona
import collections
from pathlib import Path
from typing import List
from typing import Sequence
from typing import Tuple
from typing import Union

import pandas as pd
from raven.utils import crs_sniffer, single_file_check
from shapely.geometry import shape, Point

GEO_URL = "http://boreas.ouranos.ca/geoserver/wfs"

# We store the contour of different hydrobasins domains
hybas_dir = Path(__file__).parent.parent / "data" / "hydrobasins_domains"
hybas_pat = "hybas_lake_{}_lev01_v1c.zip"

# This could be inferred from existing files in hybas_dir
hybas_regions = ["na", "ar"]

hybas_domains = {dom: hybas_dir / hybas_pat.format(dom) for dom in hybas_regions}


"""
Working assumptions for this module
-----------------------------------

* Point coordinates are passed as shapely.geometry.Point instances.
* BBox coordinates are passed as (lon1, lat1, lon2, lat2)
* Shapes (polygons) are passed as ?
* All functions that require a CRS have a CRS argument with a default set to WSG84 (?)
*

"""


def feature_contains(
    point: Tuple[Union[int, float, str], Union[str, float, int]],
    shp: Union[str, Path, List[Union[str, Path]]],
) -> Union[dict, bool]:
    """Return the first feature containing a location.

    Parameters
    ----------
    point : Tuple[Union[int, float, str], Union[str, float, int]]
      Geographic coordinates of a point (lon, lat).
    shp : Union[str, Path, List[str, Path]]
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
        point = Point(point)
    elif isinstance(point, Point):
        pass
    else:
        raise ValueError(
            "point should be shapely.Point or tuple of coordinates, got : {} of type({})".format(
                point, type(point)
            )
        )

    shape_crs = crs_sniffer(single_file_check(shp))

    if isinstance(shp, list):
        shp = shp[0]

    for i, layer_name in enumerate(fiona.listlayers(str(shp))):
        with fiona.open(shp, "r", crs=shape_crs, layer=i) as vector:
            for f in vector.filter(bbox=(point.x, point.y, point.x, point.y)):
                return f

    return False


def hydrobasins_upstream_ids(fid: str, df: pd.DataFrame) -> pd.Series:
    """Return a list of hydrobasins features located upstream.

    Parameters
    ----------
    fid : str
      HYBAS_ID of the downstream feature.
    df : pd.DataFrame
      Watershed attributes.

    Returns
    -------
    pd.Series
      Basins ids including `fid` and its upstream contributors.
    """

    def upstream_ids(bdf, bid):
        return bdf[bdf["NEXT_DOWN"] == bid]["HYBAS_ID"]

    # Locate the downstream feature
    ds = df.set_index("HYBAS_ID").loc[fid]

    # Do a first selection on the main basin ID of the downstream feature.
    sub = df[df["MAIN_BAS"] == ds["MAIN_BAS"]]

    # Find upstream basins
    up = [fid]
    for b in up:
        tmp = upstream_ids(sub, b)
        if len(tmp):
            up.extend(tmp)

    return sub[sub["HYBAS_ID"].isin(up)]


def hydrobasins_aggregate(gdf: pd.DataFrame = None) -> pd.Series:
    """Aggregate multiple hydrobasin watersheds into a single geometry.

    Parameters
    ----------
    gdf : pd.DataFrame
      Watershed attributes indexed by HYBAS_ID

    Returns
    -------
    pd.Series
    """

    # TODO: Review. Not sure it all makes sense.
    def aggfunc(x):
        if x.name in ["COAST", "DIST_MAIN", "DIST_SINK"]:
            return x.min()
        elif x.name in ["SUB_AREA", "LAKE"]:
            return x.sum()
        else:
            return x[0]

    # Buffer function to fix invalid geometries
    gdf["geometry"] = gdf.buffer(0)

    return gdf.dissolve(by="MAIN_BAS", aggfunc=aggfunc)


def get_bbox(vector: str, all_features: bool = True) -> list:
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
        with fiona.open(vector, "r") as src:
            for feature in src:
                geom = shape(feature["geometry"])
                return geom.bounds

    with fiona.open(vector, "r") as src:
        return src.bounds


def get_raster_wcs(
    coordinates: Sequence[Union[int, float, str]],
    geographic: bool = True,
    layer: str = None,
) -> bytes:
    """Return a subset of a raster image from the local GeoServer via WCS 2.0.1 protocol.

    For geoggraphic rasters, subsetting is based on WGS84 (Long, Lat) boundaries. If not geographic, subsetting based
    on projected coordinate system (Easting, Northing) boundries.

    Parameters
    ----------
    coordinates : Sequence[Union[int, float, str]]
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
        x, y = "Long", "Lat"
    else:
        x, y = "E", "N"

    wcs = WebCoverageService("http://boreas.ouranos.ca/geoserver/ows", version="2.0.1")

    try:
        resp = wcs.getCoverage(
            identifier=[layer],
            format="image/tiff",
            subsets=[(x, left, right), (y, down, up)],
        )

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


def select_hybas_domain(
    bbox: Tuple[
        Union[int, float], Union[int, float], Union[int, float], Union[int, float]
    ]
) -> str:
    """
    Provided a given coordinate or boundary box, return the domain name of the geographic region
     the coordinate is located within.

    Parameters
    ----------
    bbox : tuple
      Geographic coordinates of the bounding box (left, down, right, up).

    Returns
    -------
    str
      The domain that the coordinate falls within. Possible results: "na", "ar".
    """

    for dom, fn in hybas_domains.items():
        with open(fn, "rb") as f:
            zf = fiona.io.ZipMemoryFile(f)
            coll = zf.open(fn.stem + ".shp")
            for _ in coll.filter(bbox=bbox):
                return dom

    raise LookupError("Could not find feature containing bbox {}.".format(bbox))


def get_hydrobasins_location_wfs(
    coordinates: Tuple[
        Union[int, float, str],
        Union[str, float, int],
        Union[str, float, int],
        Union[str, float, int],
    ] = None,
    level: int = 12,
    lakes: bool = True,
    domain: str = None,
) -> str:
    """Return features from the USGS HydroBASINS data set using bounding box coordinates and WFS 1.1.0 protocol.

    For geographic rasters, subsetting is based on WGS84 (Long, Lat) boundaries. If not geographic, subsetting based
    on projected coordinate system (Easting, Northing) boundaries.

    Parameters
    ----------
    coordinates : Tuple[Union[str, float, int], Union[str, float, int], Union[str, float, int], Union[str, float, int]]
      Geographic coordinates of the bounding box (left, down, right, up).
    level : int
      Level of granularity requested for the lakes vector (1:12). Default: 12.
    lakes : bool
      Whether or not the vector should include the delimitation of lakes.
    domain : str
      The domain of the HydroBASINS data. Possible values:"na", "ar".

    Returns
    -------
    str
      A GML-encoded vector feature.

    """
    from owslib.wfs import WebFeatureService

    layer = "public:USGS_HydroBASINS_{}{}_lev{}".format(
        "lake_" if lakes else "", domain, level
    )

    if coordinates is not None:
        wfs = WebFeatureService(GEO_URL, version="1.1.0", timeout=30)
        try:
            resp = wfs.getfeature(
                typename=layer, bbox=coordinates, srsname="urn:x-ogc:def:crs:EPSG:4326"
            )
        except Exception as e:
            raise Exception(e)
    else:
        raise NotImplementedError

    data = resp.read()
    return data


# TODO: Fix docstring. This does not return features, but a URL which will trigger a remote service returning features.
def get_hydrobasins_attributes_wfs(
    attribute: str = None,
    value: Union[str, float, int] = None,
    level: int = 12,
    lakes: bool = True,
    domain: str = None,
) -> str:
    """Return features from the USGS HydroBASINS data set using attribute value selection and WFS 1.1.0 protocol.

    For geographic rasters, subsetting is based on WGS84 (Long, Lat) boundaries. If not geographic, subsetting based
    on projected coordinate system (Easting, Northing) boundaries.

    Parameters
    ----------
    attribute : str
      Attribute/field to be queried.
    value: Union[str, float, int]
      Value for attribute queried.
    level : int
      Level of granularity requested for the lakes vector (1:12). Default: 12.
    lakes : bool
      Whether or not the vector should include the delimitation of lakes.
    domain : str
      The domain of teh HydroBASINS data. Possible values:"na", "ar".

    Returns
    -------
    str
      URL to the GeoJSON-encoded WFS response.

    """
    from requests import Request
    from owslib.fes import PropertyIsLike
    from lxml import etree

    layer = "public:USGS_HydroBASINS_{}{}_lev{}".format(
        "lake_" if lakes else "", domain, level
    )

    if attribute is not None and value is not None:

        try:
            attribute = str(attribute)
            value = str(value)

        except ValueError:
            raise Exception("Unable to cast attribute/filter to string")

        try:
            filter_request = PropertyIsLike(
                propertyname=attribute, literal=value, wildCard="*"
            )
            filterxml = etree.tostring(filter_request.toXML()).decode("utf-8")
            params = dict(
                service="WFS",
                version="1.1.0",
                request="GetFeature",
                typename=layer,
                outputFormat="json",
                filter=filterxml,
            )

            q = Request("GET", GEO_URL, params=params).prepare().url

        except Exception as e:
            raise Exception(e)

    else:
        raise NotImplementedError

    return q
