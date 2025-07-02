import logging
import tempfile
from pathlib import Path
from typing import Optional, Union

import numpy as np
import rasterio
from osgeo.gdal import Dataset, DEMProcessing
from shapely.geometry import GeometryCollection, MultiPolygon, Polygon, shape

from raven.utilities.geo import generic_raster_clip

# See: https://kokoalberti.com/articles/geotiff-compression-optimization-guide/
# or 'compress=deflate' or 'compress=zstd' or 'compress=lerc' or others
GDAL_TIFF_COMPRESSION_OPTION = "compress=lzw"

LOGGER = logging.getLogger("RavenPy")


def geom_prop(geom: Union[Polygon, MultiPolygon, GeometryCollection]) -> dict:
    """
    Return a dictionary of geometry properties.

    Parameters
    ----------
    geom : Polygon or MultiPolygon or GeometryCollection
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
    gravelius = length / 2 / np.sqrt(np.pi * area)

    # JSON serializable output
    properties = {
        "area": area,
        "centroid": (lon, lat),
        "perimeter": length,
        "gravelius": float(gravelius),
    }
    return properties


def dem_prop(
    dem: Union[str, Path],
    geom: Union[Polygon, MultiPolygon, list[Union[Polygon, MultiPolygon]]] = None,
    directory: Union[str, Path] = None,
) -> dict[str, float]:
    """Return raster properties for each geometry.

    This

    Parameters
    ----------
    dem : str or Path
        DEM raster in reprojected coordinates.
    geom : Polygon or MultiPolygon or List[Polygon or MultiPolygon]
        Geometry over which aggregate properties will be computed. If None compute properties over entire raster.
    directory : str or Path
        Folder to save the GDAL terrain analysis outputs.

    Returns
    -------
    dict
        Dictionary storing mean elevation [m], slope [deg] and aspect [deg] as float.
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
    slope = gdal_slope_analysis(fns["dem"], set_output=fns["slope"])

    # Compute aspect
    aspect = gdal_aspect_analysis(fns["dem"], set_output=fns["aspect"])
    aspect_mean = circular_mean_aspect(aspect)

    # JSON serializable output
    properties = {
        "elevation": float(elevation.mean()),
        "slope": float(slope.mean()),
        "aspect": float(aspect_mean),
    }
    return properties


def gdal_slope_analysis(
    dem: Union[str, Path],
    set_output: Optional[Union[str, Path]] = None,
    units: str = "degree",
) -> np.ndarray:
    """Return the slope of the terrain from the DEM.

    The slope is the magnitude of the gradient of the elevation.

    Parameters
    ----------
    dem : str or Path
        Path to file storing DEM.
    set_output : str or Path, optional
        If set to a valid filepath, will write to this path, otherwise will use an in-memory gdal.Dataset.
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
    if isinstance(dem, Path):
        dem = str(dem)
    if set_output:
        if isinstance(set_output, (str, Path)):
            set_output = str(set_output)
            DEMProcessing(
                set_output,
                dem,
                "slope",
                slopeFormat=units,
                format="GTiff",
                band=1,
                creationOptions=[GDAL_TIFF_COMPRESSION_OPTION],
            )
            with rasterio.open(set_output) as src:
                return np.ma.masked_values(src.read(1), value=-9999)
        else:
            raise ValueError()
    else:
        set_output = DEMProcessing(
            "",
            dem,
            "slope",
            slopeFormat=units,
            format="MEM",
            band=1,
        )
        return np.ma.masked_values(set_output.ReadAsArray(), value=-9999)


def gdal_aspect_analysis(
    dem: Union[str, Path],
    set_output: Union[str, Path, bool] = False,
    flat_values_are_zero: bool = False,
) -> Union[np.ndarray, Dataset]:
    """Return the aspect of the terrain from the DEM.

    The aspect is the compass direction of the steepest slope (0: North, 90: East, 180: South, 270: West).

    Parameters
    ----------
    dem : str or Path
        Path to file storing DEM.
    set_output : str or Path or bool
        If set to a valid filepath, will write to this path, otherwise will use an in-memory gdal.Dataset.
    flat_values_are_zero : bool
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
    if isinstance(dem, Path):
        dem = str(dem)
    if set_output:
        if isinstance(set_output, (str, Path)):
            set_output = str(set_output)
            DEMProcessing(
                destName=set_output,
                srcDS=dem,
                processing="aspect",
                zeroForFlat=flat_values_are_zero,
                format="GTiff",
                band=1,
                creationOptions=[GDAL_TIFF_COMPRESSION_OPTION],
            )
            with rasterio.open(set_output) as src:
                return np.ma.masked_values(src.read(1), value=-9999)
        else:
            raise ValueError()

    else:
        set_output = DEMProcessing(
            destName="",
            srcDS=dem,
            processing="aspect",
            zeroForFlat=flat_values_are_zero,
            format="MEM",
            band=1,
        )
        return np.ma.masked_values(set_output.ReadAsArray(), value=-9999)


def circular_mean_aspect(angles: np.ndarray) -> np.ndarray:
    """
    Return the mean angular aspect based on circular arithmetic approach.

    Parameters
    ----------
    angles : np.ndarray
        Array of aspect angles.

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
