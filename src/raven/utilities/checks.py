"""Checks for various geospatial and IO conditions."""

import collections
import logging
import warnings
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import fiona
import rasterio
from pyproj import CRS
from pyproj.exceptions import CRSError
from raven.utilities import io
from shapely.geometry import GeometryCollection, MultiPolygon, Point, shape
from shapely.geometry.base import BaseGeometry


logger = logging.getLogger(__name__)


def single_file_check(file_list: Sequence[str | Path]) -> Any:
    """
    Return the first element of a file list and raise an error if the list is empty or contains more than one element.

    Parameters
    ----------
    file_list : Sequence of str or Path
        A list of files.
    """
    if isinstance(file_list, (str, Path)):
        return file_list

    try:
        if len(file_list) > 1:
            msg = "Multi-file handling for file is not supported. Exiting."
            raise NotImplementedError(msg)
        elif len(file_list) == 0:
            msg = "No files found. Exiting."
            raise FileNotFoundError(msg)
        return file_list[0]
    except (FileNotFoundError, NotImplementedError) as e:
        logger.error(e)
        raise


def boundary_check(
    *args: str | Path,
    max_y: float = 60,
    min_y: float = -60,
) -> None:
    r"""
    Verify that boundaries do not exceed specific latitudes for geographic coordinate data.

    Parameters
    ----------
    *args : Sequence of str or Path
        str or Path to files.
    max_y : int or float
        Maximum value allowed for latitude. Default: 60.
    min_y : int or float
        Minimum value allowed for latitude. Default: -60.

    Warnings
    --------
    UserWarning
        If boundaries exceed specified latitude boundaries.
    """
    vectors = (".gml", ".shp", ".geojson", ".gpkg", ".json")
    rasters = (".tif", ".tiff")

    for file in args:
        file = Path(file)

        try:
            if file.suffix.lower() in vectors:
                src = fiona.open(file, "r")
            elif file.suffix.lower() in rasters:
                src = rasterio.open(file, "r")
            else:
                raise FileNotFoundError()

            try:
                geographic = CRS(src.crs).is_geographic
            except CRSError:
                geographic = True
            src_min_y, src_max_y = src.bounds[1], src.bounds[3]
            if geographic and (src_max_y > max_y or src_min_y < min_y):
                msg = f"Vector {file} contains geometries in high latitudes. Verify choice of projected CRS is appropriate for analysis."
                logger.warning(msg)
                warnings.warn(msg, UserWarning)
            if not geographic:
                msg = f"Vector {file} is not in a geographic coordinate system."
                logger.warning(msg)
                warnings.warn(msg, UserWarning)
            src.close()

        except FileNotFoundError:
            msg = f"Unable to read boundaries from {file}"
            logger.error(msg)
            raise


def multipolygon_check(geom: GeometryCollection) -> None:
    """
    Perform a check to verify a geometry is a MultiPolygon.

    Parameters
    ----------
    geom : GeometryCollection
        The geometry to check.
    """
    if not isinstance(geom, GeometryCollection):
        try:
            geom = shape(geom)
        except AttributeError:
            logger.error("Unable to load argument as shapely.geometry.shape().")
            raise

    if isinstance(geom, MultiPolygon):
        logger.warning("Shape is a Multipolygon.")


PointType = tuple[int | float | str, int | float | str] | BaseGeometry


def feature_contains(
    point: PointType,
    shp: str | Path | list[str | Path],
) -> dict | bool:
    """
    Return the first feature containing a location.

    Parameters
    ----------
    point : tuple of [int or float or str, int or float or str] or Point
        Geographic coordinates of a point (lon, lat) or a shapely Point.
    shp : str or Path or list of str or Path
        The path to the file storing the geometries.

    Returns
    -------
    dict or bool
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
        raise TypeError(
            f"point should be shapely.Point or tuple of coordinates, got : {point} of type({type(point)})"
        )

    shape_crs = io.crs_sniffer(single_file_check(shp))

    if isinstance(shp, list):
        shp = shp[0]

    for i, layer_name in enumerate(fiona.listlayers(str(shp))):
        with fiona.open(shp, "r", crs=shape_crs, layer=i) as vector:
            for f in vector.filter(bbox=(point.x, point.y, point.x, point.y)):
                return dict(f)

    return False
