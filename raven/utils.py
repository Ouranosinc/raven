import json
import logging
import math
import re
import tarfile
import tempfile
import zipfile
from functools import partial
from pathlib import Path
from re import search
from typing import Any
from typing import List
from typing import Optional
from typing import Sequence
from typing import Tuple
from typing import Union

import fiona
import fiona.crs
import numpy as np
import pyproj
import rasterio
import rasterio.mask
import rasterio.vrt
import rasterio.warp
from gdal import DEMProcessing
from rasterio.crs import CRS
from shapely.geometry import GeometryCollection
from shapely.geometry import MultiPolygon
from shapely.geometry import Polygon
from shapely.geometry import mapping
from shapely.geometry import shape
from shapely.ops import transform

LOGGER = logging.getLogger("RAVEN")

# See: https://kokoalberti.com/articles/geotiff-compression-optimization-guide/
GDAL_TIFF_COMPRESSION_OPTION = "compress=lzw"  # or 'compress=deflate' or 'compress=zstd' or 'compress=lerc' or others
RASTERIO_TIFF_COMPRESSION = "lzw"

WGS84 = "+init=epsg:4326"
WGS84_PROJ4 = "+proj=longlat +datum=WGS84 +no_defs"
LAEA = "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
WORLDMOLL = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
ALBERS_NAM = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"


def address_append(address: Union[str, Path]) -> str:
    """
    Formats a URL/URI to be more easily read with libraries such as "rasterstats"

    Parameters
    ----------
    address: Union[str, Path]
      URL/URI to a potential zip or tar file

    Returns
    -------
    str
      URL/URI prefixed for archive type
    """
    zipped = search(r"(\.zip)", str(address))
    tarred = search(r"(\.tar)", str(address))

    try:
        if zipped:
            return "zip://{}".format(address)
        elif tarred:
            return "tar://{}".format(address)
        else:
            LOGGER.info("No changes made to address.")
            return address
    except Exception as e:
        msg = "Failed to prefix or parse URL {}: {}".format(address, e)
        LOGGER.error(msg)


def generic_extract_archive(
    resources: Union[str, Path, List[Union[bytes, str, Path]]],
    output_dir: Optional[Union[str, Path]] = None,
) -> List[str]:
    """Extract archives (tar/zip) to a working directory.

    Parameters
    ----------
    resources: Union[str, Path, List[Union[bytes, str, Path]]]
      list of archive files (if netCDF files are in list, they are passed and returned as well in the return).
    output_dir: Optional[Union[str, Path]]
      string or Path to a working location (default: temporary folder).

    Returns
    -------
    list
      List of original or of extracted files
    """

    archive_types = [".tar", ".zip", ".7z"]
    output_dir = output_dir or tempfile.gettempdir()

    if not isinstance(resources, list):
        resources = [resources]

    files = []

    for arch in resources:
        if any(ext in str(arch).lower() for ext in archive_types):
            try:
                LOGGER.debug("archive=%s", arch)
                file = Path(arch).name

                if file.endswith(".nc"):
                    files.append(Path(output_dir.join(arch)))
                elif file.endswith(".tar"):
                    with tarfile.open(arch, mode="r") as tar:
                        tar.extractall(path=output_dir)
                        files.extend(
                            [str(Path(output_dir).joinpath(f)) for f in tar.getnames()]
                        )
                elif file.endswith(".zip"):
                    with zipfile.ZipFile(arch, mode="r") as zf:
                        zf.extractall(path=output_dir)
                        files.extend(
                            [str(Path(output_dir).joinpath(f)) for f in zf.namelist()]
                        )
                elif file.endswith(".7z"):
                    msg = "7z file extraction is not supported at this time"
                    LOGGER.warning(msg)
                    raise UserWarning(msg)
                else:
                    LOGGER.debug('File extension "{}" unknown'.format(file))
            except Exception as e:
                LOGGER.error("Failed to extract sub archive {}: {}".format(arch, e))
        else:
            LOGGER.warning("No archives found. Continuing...")
            return resources

    return files


def archive_sniffer(
    archives: Union[str, Path, List[Union[str, Path]]],
    working_dir: Union[str, Path],
    extensions: Optional[List[str]] = None,
) -> List[Union[str, Path]]:
    """Return a list of locally unarchived files that match the desired extensions.

    Parameters
    ----------
    archives : Union[str, Path, List[Union[str, Path]]]
      archive location or list of archive locations
    working_dir : Union[str, path]
      string or Path to a working location
    extensions : List[str]
      list of accepted extensions

    Returns
    -------
    List[Union[str, Path]]
      List of files with matching accepted extensions
    """
    potential_files = []

    if not extensions:
        extensions = [".gml", ".shp", ".geojson", ".gpkg", ".json"]

    decompressed_files = generic_extract_archive(archives, output_dir=working_dir)
    for file in decompressed_files:
        if any(ext in Path(file).suffix for ext in extensions):
            potential_files.append(file)
    return potential_files


def crs_sniffer(*args: Sequence[Union[str, Path]]) -> Union[List[str], str]:
    """Return the list of CRS found in files.

    Parameters
    ----------
    args : Sequence[Union[str, Path]]
      Paths to the files to examine.

    Returns
    -------
    Union[List[str], str]
      Returns either a list of CRSes or a single CRS definition, depending on the number of instances found.
    """
    crs_list = []
    vectors = (".gml", ".shp", ".geojson", ".gpkg", ".json")
    rasters = (".tif", ".tiff")

    for file in args:
        found_crs = False
        try:
            if str(file).lower().endswith(vectors):
                if str(file).lower().endswith(".gpkg"):
                    if len(fiona.listlayers(file)) > 1:
                        raise NotImplementedError
                with fiona.open(file, "r") as src:
                    found_crs = fiona.crs.to_string(src.crs)
            elif str(file).lower().endswith(rasters):
                with rasterio.open(file, "r") as src:
                    found_crs = CRS(src.crs).to_proj4()
            else:
                raise FileNotFoundError("Invalid filename suffix")
        except FileNotFoundError as e:
            msg = "{}: Unable to open file {}".format(e, args)
            LOGGER.warning(msg)
            raise Exception(msg)
        except NotImplementedError as e:
            msg = "{}: Multilayer GeoPackages are currently unsupported".format(e)
            LOGGER.error(msg)
            raise Exception(msg)
        except RuntimeError:
            pass

        crs_list.append(found_crs)

    if crs_list is None:
        msg = "No CRS definitions found in {}.".format(args)
        raise FileNotFoundError(msg)

    if len(crs_list) == 1:
        if crs_list[0] == "":
            msg = "No CRS definitions found in {}. Using {}".format(args, WGS84_PROJ4)
            LOGGER.warning(msg)
            return WGS84_PROJ4
        return crs_list[0]
    return crs_list


def raster_datatype_sniffer(file: Union[str, Path]) -> str:
    """Return the type of the raster stored in the file.

    Parameters
    ----------
    file : Union[str, Path]
      Path to file.

    Returns
    -------
    str
      rasterio datatype of array values
    """
    try:
        with rasterio.open(file, "r") as src:
            dtype = src.dtypes[0]
        return dtype
    except Exception as e:
        msg = "{}: Unable to read data type from {}".format(e, file)
        LOGGER.exception(msg)
        raise ValueError(msg)


def parse_lonlat(lonlat: Union[str, Tuple[str, str]]) -> Tuple[float, float]:
    """Return longitude and latitude from a string.

    Parameters
    ----------
    lonlat : Union[str, Tuple[str, str]]
      A tuple or a str of lon and lat coordinates.

    Returns
    -------
    Tuple[float, float]
    """
    try:
        if isinstance(lonlat, str):
            lon, lat = tuple(map(float, re.findall(r"[-+]?[0-9]*\.?[0-9]+", lonlat)))
        elif isinstance(lonlat, tuple):
            lon, lat = map(float, lonlat)
        else:
            raise ValueError
        return lon, lat
    except Exception as e:
        msg = "Failed to parse longitude, latitude coordinates {}".format(lonlat)
        raise Exception(msg) from e


def single_file_check(file_list: List[Union[str, Path]]) -> Any:
    """Return the first element of a file list. Raise an error if the list is empty or contains more than one element.

    Parameters
    ----------
    file_list : List[Union[str, Path]]
    """
    try:
        if len(file_list) > 1:
            msg = "Multi-file handling for file is not supported. Exiting."
            raise NotImplementedError(msg)
        elif len(file_list) == 0:
            msg = "No files found. Exiting."
            raise FileNotFoundError(msg)
        return file_list[0]
    except Exception as e:
        msg = "{}: Unspecified error. Exiting,"
        LOGGER.error(msg.format(e))
        raise Exception(msg)


def boundary_check(
    *args: Sequence[Union[str, Path]],
    max_y: Union[int, float] = 60,
    min_y: Union[int, float] = -60
) -> None:
    """Verify that boundaries do not exceed specific latitudes for geographic coordinate data. Raise a warning if so.

    Parameters
    ----------
    *args : Sequence[Union[str, Path]]
      listing of strings or paths to files
    max_y : Union[int, float]
      Maximum value allowed for latitude. Default: 60.
    min_y : Union[int, float]
      Minimum value allowed for latitude. Default: -60.
    """
    vectors = (".gml", ".shp", ".geojson", ".gpkg", ".json")
    rasters = (".tif", ".tiff")
    for file in args:
        src = None
        try:
            if str(file).lower().endswith(vectors):
                src = fiona.open(file, "r")
            elif str(file).lower().endswith(rasters):
                src = rasterio.open(file, "r")
            else:
                FileNotFoundError("Invalid filename suffix")

            geographic = CRS(src.crs).is_geographic
            if geographic and (src.bounds > max_y or src.bounds < min_y):
                msg = (
                    "Vector {} contains geometries in high latitudes."
                    " Verify choice of projected CRS is appropriate for analysis.".format(
                        file
                    )
                )
                LOGGER.warning(msg)
                UserWarning(msg)
            src.close()

        except Exception as e:
            msg = "{}: Unable to read boundaries from {}".format(e, args)
            LOGGER.exception(msg)
    return


def multipolygon_check(geom: GeometryCollection) -> None:
    """Perform a check to verify a geometry is a MultiPolygon

    Parameters
    ----------
    geom : GeometryCollection

    Returns
    -------
    None
    """
    if not isinstance(type(geom), GeometryCollection):
        try:
            geom = shape(geom)
        except Exception as e:
            LOGGER.error(
                "{}: Unable to load vector as shapely.geometry.shape().".format(e)
            )

    if isinstance(type(geom), MultiPolygon):
        LOGGER.warning("Shape is a Multipolygon.")
    return


def geom_transform(
    geom: GeometryCollection,
    source_crs: Union[str, CRS] = WGS84,
    target_crs: Union[str, CRS] = None,
) -> GeometryCollection:
    """Change the projection of a geometry.

    Assuming a geometry's coordinates are in a `source_crs`, compute the new coordinates under the `target_crs`.

    Parameters
    ----------
    geom : GeometryCollection
      Source geometry.
    source_crs : Union[str, CRS]
      Projection identifier (proj4) for the source geometry, e.g. '+proj=longlat +datum=WGS84 +no_defs'.
    target_crs : Union[str, CRS]
      Projection identifier (proj4) for the target geometry.

    Returns
    -------
    GeometryCollection
      Reprojected geometry.
    """
    try:
        reprojected = transform(
            partial(pyproj.transform, pyproj.Proj(source_crs), pyproj.Proj(target_crs)),
            geom,
        )
        return reprojected
    except Exception as e:
        msg = "{}: Failed to reproject geometry".format(e)
        LOGGER.error(msg)
        raise Exception(msg)


def geom_prop(geom: Union[Polygon, MultiPolygon, GeometryCollection]) -> dict:
    """Return a dictionary of geometry properties.

    Parameters
    ----------
    geom : Union[Polygon, MultiPolygon, GeometryCollection]
      Geometry to analyze.

    Returns
    -------
    dict
      Dictionary storing polygon area, centroid location, perimeter and gravelius shape index.

    Notes
    -----
    Some of the properties should be computed using an equal-area projection.
    """

    geom = shape(geom)
    lon, lat = geom.centroid.x, geom.centroid.y
    if (lon > 180) or (lon < -180) or (lat > 90) or (lat < -90):
        LOGGER.warning("Shape centroid is not in decimal degrees.")
    area = geom.area
    length = geom.length
    gravelius = length / 2 / math.sqrt(math.pi * area)
    parameters = {
        "area": area,
        "centroid": (lon, lat),
        "perimeter": length,
        "gravelius": gravelius,
    }
    return parameters


def dem_prop(
    dem: Union[str, Path],
    geom: Union[Polygon, MultiPolygon, List[Union[Polygon, MultiPolygon]]] = None,
    directory: Union[str, Path] = None,
) -> dict:
    """Return raster properties for each geometry.

    This

    Parameters
    ----------
    dem : Union[str, Path]
      DEM raster in reprojected coordinates.
    geom : Union[Polygon, MultiPolygon, List[Union[Polygon, MultiPolygon]]]
      Geometry over which aggregate properties will be computed. If None compute properties over entire raster.
    directory : Union[str, Path]
      Folder to save the GDAL terrain anlaysis outputs

    Returns
    -------
    dict
      Dictionary storing mean elevation [m], slope [deg] and aspect [deg].
    """

    fns = dict()
    fns["dem"] = (
        tempfile.NamedTemporaryFile(
            prefix="dem", suffix=".tiff", dir=directory, delete=False
        ).name
        if geom is not None
        else dem
    )
    for key in ["slope", "aspect"]:
        fns[key] = tempfile.NamedTemporaryFile(
            prefix=key, suffix=".tiff", dir=directory, delete=False
        ).name

    # Clip to relevant area or read original raster
    if geom is None:
        with rasterio.open(dem) as f:
            elevation = f.read(1, masked=True)
    else:
        generic_raster_clip(raster=dem, output=fns["dem"], geometry=geom)
        with rasterio.open(fns["dem"]) as f:
            elevation = f.read(1, masked=True)

    # Compute slope
    slope = gdal_slope_analysis(fns["dem"], output=fns["slope"])

    # Compute aspect
    aspect = gdal_aspect_analysis(fns["dem"], output=fns["aspect"])
    aspect_mean = circular_mean_aspect(aspect)

    return {"elevation": elevation.mean(), "slope": slope.mean(), "aspect": aspect_mean}


# It's a bit weird to have to pass the output file name as an argument, since you return an in-memory array.
# Can you keep everything in memory ?

# Response: The DEMPROCESSING function is basically a thin GDAL wrapper for a os.subprocess call. There is
# unfortunately no way to call the gdal slope/aspect calculation from Python directly so it demands an output filename.
# Since it technically is writing the information to the file, this function could use a generic named temporary file
# to perform the analysis and return the array all in memory I suppose. I've added this as an option here.


def gdal_slope_analysis(
    dem: Union[str, Path], output: Union[str, Path] = None, units: str = "degree"
) -> np.ndarray:
    """Return the slope of the terrain from the DEM.

    The slope is the magnitude of the gradient of the elevation.

    Parameters
    ----------
    dem : Union[str, Path]
      Path to file storing DEM.
    output : Union[str, Path]
      Path to output file.
    units : str
      Slope units. Default: 'degree'.

    Returns
    -------
    np.ndarray
      Slope array.

    Notes
    -----
    Ensure that the DEM is in a *projected coordinate*, not a geographic coordinate system, so that the
    horizontal scale is the same as the vertical scale (m).

    """
    if output is None:
        output = tempfile.NamedTemporaryFile().name
    if isinstance(dem, Path):
        dem = str(dem)
    if isinstance(output, Path):
        output = str(output)

    DEMProcessing(
        output,
        dem,
        "slope",
        slopeFormat=units,
        format="GTiff",
        band=1,
        creationOptions=[GDAL_TIFF_COMPRESSION_OPTION],
    )
    with rasterio.open(output) as src:
        return np.ma.masked_values(src.read(1), value=-9999)


def gdal_aspect_analysis(
    dem: Union[str, Path],
    output: Union[str, Path] = None,
    flat_values_are_zero: bool = False,
) -> np.ndarray:
    """Return the aspect of the terrain from the DEM.

    The aspect is the compass direction of the steepest slope (0: North, 90: East, 180: South, 270: West).

    Parameters
    ----------
    dem : Union[str, Path]
      Path to file storing DEM.
    output : Union[str, Path]
      Path to output file.
    flat_values_are_zero: bool
      Designate flat values with value zero. Default: -9999.

    Returns
    -------
    np.ndarray
      Aspect array.

    Notes
    -----
    Ensure that the DEM is in a *projected coordinate*, not a geographic coordinate system, so that the
    horizontal scale is the same as the vertical scale (m).
    """
    if output is None:
        output = tempfile.NamedTemporaryFile().name
    if isinstance(dem, Path):
        dem = str(dem)
    if isinstance(output, Path):
        output = str(output)

    DEMProcessing(
        destName=output,
        srcDS=dem,
        processing="aspect",
        zeroForFlat=flat_values_are_zero,
        format="GTiff",
        band=1,
        creationOptions=[GDAL_TIFF_COMPRESSION_OPTION],
    )
    with rasterio.open(output) as src:
        return np.ma.masked_values(src.read(1), value=-9999)


def circular_mean_aspect(angles: np.ndarray) -> np.ndarray:
    """Return the mean angular aspect based on circular arithmetic approach

    Parameters
    ----------
    angles: np.ndarray
      Array of aspect angles

    Returns
    -------
    np.ndarray
      Circular mean of aspect array.
    """
    # Circular statistics needed for mean angular aspect
    # Example from: https://gis.stackexchange.com/a/147135/65343

    n = len(angles)
    sine_mean = np.divide(np.sum(np.sin(np.radians(np.ma.masked_array(angles)))), n)
    cosine_mean = np.divide(np.sum(np.cos(np.radians(np.ma.masked_array(angles)))), n)
    vector_mean = np.arctan2(sine_mean, cosine_mean)
    degrees = np.degrees(vector_mean)

    if degrees < 0:
        return degrees + 360
    return degrees


def generic_raster_clip(
    raster: Union[str, Path],
    output: Union[str, Path],
    geometry: Union[Polygon, MultiPolygon, List[Union[Polygon, MultiPolygon]]],
    touches: bool = False,
    fill_with_nodata: bool = True,
    padded: bool = True,
    raster_compression: str = RASTERIO_TIFF_COMPRESSION,
) -> None:
    """
    Crop a raster file to a given geometry.

    Parameters
    ----------
    raster : Union[str, Path]
      Path to input raster.
    output : Union[str, Path]
      Path to output raster.
    geometry : Union[Polygon, MultiPolygon, List[Union[Polygon, MultiPolygon]]
      Geometry defining the region to crop.
    touches : bool
      Whether or not to include cells that intersect the geometry. Default: True.
    fill_with_nodata: bool
      Whether or not to keep pixel values for regions outside of shape or set as nodata. Default: True.
    padded: bool
      Whether or not to add a half-pixel buffer to shape before masking
    raster_compression : str
      Level of data compression. Default: 'lzw'.

    Returns
    -------
    None
    """
    if not (type(geometry) in (list, tuple)):
        geometry = [geometry]

    with rasterio.open(raster, "r") as src:
        mask_image, mask_affine = rasterio.mask.mask(
            src,
            geometry,
            crop=True,
            pad=padded,
            all_touched=touches,
            filled=fill_with_nodata,
        )
        mask_meta = src.meta.copy()
        mask_meta.update(
            {
                "driver": "GTiff",
                "height": mask_image.shape[1],
                "width": mask_image.shape[2],
                "transform": mask_affine,
                "compress": raster_compression,
            }
        )

        # Write the new masked image
        with rasterio.open(output, "w", **mask_meta) as dst:
            dst.write(mask_image)
    return


def generic_raster_warp(
    raster: Union[str, Path],
    output: Union[str, Path],
    target_crs: Union[str, dict, CRS],
    raster_compression: str = RASTERIO_TIFF_COMPRESSION,
) -> None:
    """
    Reproject a raster file.

    Parameters
    ----------
    raster : Union[str, Path]
      Path to input raster.
    output : Union[str, Path]
      Path to output raster.
    target_crs : str or dict
      Target projection identifier.
    raster_compression: str
      Level of data compression. Default: 'lzw'.

    Returns
    -------
    None
    """
    with rasterio.open(raster, "r") as src:
        # Reproject raster using WarpedVRT class
        with rasterio.vrt.WarpedVRT(src, crs=target_crs) as vrt:
            # Calculate grid properties based on projection
            affine, width, height = rasterio.warp.calculate_default_transform(
                src.crs, target_crs, src.width, src.height, *src.bounds
            )

            # Copy relevant metadata from parent raster
            metadata = src.meta.copy()
            metadata.update(
                {
                    "driver": "GTiff",
                    "height": height,
                    "width": width,
                    "transform": affine,
                    "crs": target_crs,
                    "compress": raster_compression,
                }
            )
            data = vrt.read()

            with rasterio.open(output, "w", **metadata) as dst:
                dst.write(data)
    return


def generic_vector_reproject(
    vector: Union[str, Path],
    projected: Union[str, Path],
    source_crs: Union[str, CRS] = WGS84_PROJ4,
    target_crs: Union[str, CRS] = None,
) -> None:
    """Reproject all features and layers within a vector file and return a GeoJSON

    Parameters
    ----------
    vector : Union[str, Path]
      Path to a file containing a valid vector layer.
    projected: Union[str, Path]
      Path to a file to be written.
    source_crs : Union[str, dict, CRS]
      Projection identifier (proj4) for the source geometry, Default: '+proj=longlat +datum=WGS84 +no_defs'.
    target_crs : Union[str, dict, CRS]
      Projection identifier (proj4) for the target geometry.

    Returns
    -------
    None
    """

    if target_crs is None:
        msg = "No target CRS is defined."
        raise ValueError(msg)

    output = {"type": "FeatureCollection", "features": []}

    if isinstance(vector, Path):
        vector = str(vector)

    for i, layer_name in enumerate(fiona.listlayers(vector)):
        with fiona.open(vector, "r", layer=i) as src:
            with open(projected, "w") as sink:
                for feature in src:
                    # Perform vector reprojection using Shapely on each feature
                    try:
                        geom = shape(feature["geometry"])
                        transformed = geom_transform(geom, source_crs, target_crs)
                        feature["geometry"] = mapping(transformed)
                        output["features"].append(feature)
                    except Exception as e:
                        msg = "{}: Unable to reproject feature {}".format(e, feature)
                        LOGGER.exception(msg)

                sink.write("{}".format(json.dumps(output)))
    return
