"""Tools for performing geospatial translations and transformations."""

import collections
import logging
from pathlib import Path
from typing import Optional, Union

import fiona
import geopandas
import pandas as pd
import rasterio
import rasterio.mask
import rasterio.vrt
import rasterio.warp
from fiona import Feature
from pyproj import CRS
from shapely.geometry import (
    GeometryCollection,
    MultiPolygon,
    Point,
    Polygon,
    mapping,
    shape,
)
from shapely.ops import transform

RASTERIO_TIFF_COMPRESSION = "lzw"
LOGGER = logging.getLogger("RavenPy")
WGS84 = 4326


def geom_transform(
    geom: Union[GeometryCollection, shape],
    source_crs: Union[str, int, CRS] = WGS84,
    target_crs: Union[str, int, CRS] = None,
) -> GeometryCollection:
    """
    Change the projection of a geometry.

    Assuming a geometry's coordinates are in a `source_crs`, compute the new coordinates under the `target_crs`.

    Parameters
    ----------
    geom : Union[GeometryCollection, shape]
        Source geometry.
    source_crs : Union[str, int, CRS]
        Projection identifier (proj4) for the source geometry, e.g. '+proj=longlat +datum=WGS84 +no_defs'.
    target_crs : Union[str, int, CRS]
        Projection identifier (proj4) for the target geometry.

    Returns
    -------
    GeometryCollection
        Reprojected geometry.
    """
    try:
        from pyproj import Transformer

        if isinstance(source_crs, int or str):
            source = CRS.from_epsg(source_crs)
        else:
            source = source_crs

        if isinstance(target_crs, int or str):
            target = CRS.from_epsg(target_crs)
        else:
            target = target_crs

        transform_func = Transformer.from_crs(source, target, always_xy=True)
        reprojected = transform(transform_func.transform, geom)

        return reprojected
    except Exception as err:
        msg = f"{err}: Failed to reproject geometry"
        LOGGER.error(msg)
        raise Exception(msg)


def generic_raster_clip(
    raster: Union[str, Path],
    output: Union[str, Path],
    geometry: Union[Polygon, MultiPolygon, list[Union[Polygon, MultiPolygon]]],
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
        Whether to include cells that intersect the geometry or not. Default: True.
    fill_with_nodata : bool
        Whether to keep pixel values for regions outside shape or set as nodata or not. Default: True.
    padded : bool
        Whether to add a half-pixel buffer to shape before masking or not. Default: True.
    raster_compression : str
        Level of data compression. Default: 'lzw'.
    """
    if not isinstance(geometry, collections.abc.Iterable):
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
    raster_compression : str
        Level of data compression. Default: 'lzw'.
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


def generic_vector_reproject(
    vector: Union[str, Path],
    projected: Union[str, Path],
    source_crs: Union[str, CRS] = WGS84,
    target_crs: Union[str, CRS] = None,
) -> None:
    """
    Reproject all features and layers within a vector file and return a GeoJSON.

    Parameters
    ----------
    vector : Union[str, Path]
        Path to a file containing a valid vector layer.
    projected : Union[str, Path]
        Path to a file to be written.
    source_crs : Union[str, dict, CRS]
        Projection identifier (proj4) for the source geometry, Default: '+proj=longlat +datum=WGS84 +no_defs'.
    target_crs : Union[str, dict, CRS]
        Projection identifier (proj4) for the target geometry.
    """

    if target_crs is None:
        raise ValueError("No target CRS is defined.")

    if isinstance(vector, Path):
        vector = vector.as_posix()

    for i, layer_name in enumerate(fiona.listlayers(vector)):
        with fiona.open(vector, "r", layer=i) as src:
            # Increase the size allocation of integer fields
            for entry, formatting in src.schema["properties"].items():
                if formatting.startswith("int"):
                    src.schema["properties"][entry] = "int64"
            with fiona.open(
                projected, "w", driver="GeoJSON", schema=src.schema, crs=target_crs
            ) as sink:
                for feature in src:
                    # Perform vector reprojection using Shapely on each feature
                    try:
                        geom = shape(feature["geometry"])
                        transformed = geom_transform(geom, source_crs, target_crs)
                        new_geom = fiona.Geometry().from_dict(mapping(transformed))
                        sink.write(
                            Feature(geometry=new_geom, properties=feature.properties)
                        )
                    except Exception as err:
                        LOGGER.exception(
                            f"{err}: Unable to reproject feature {feature}"
                        )
                        raise


def determine_upstream_ids(
    fid: str,
    df: Union[pd.DataFrame, geopandas.GeoDataFrame],
    basin_field: str = None,
    downstream_field: str = None,
    basin_family: Optional[str] = None,
) -> Union[pd.DataFrame, geopandas.GeoDataFrame]:
    """
    Return a list of upstream features by evaluating the downstream networks.

    Parameters
    ----------
    fid : str
        feature ID of the downstream feature of interest.
    df : pd.DataFrame
        A Dataframe comprising the watershed attributes.
    basin_field : str
        The field used to determine the id of the basin according to hydro project.
    downstream_field : str
        The field identifying the downstream subbasin for the hydro project.
    basin_family : str, optional
        Regional watershed code (For HydroBASINS dataset).

    Returns
    -------
    pandas.DataFrame or geopandas.GeoDataFrame
        Basins ids including `fid` and its upstream contributors.
    """

    def upstream_ids(bdf, bid):
        return bdf[bdf[downstream_field] == bid][basin_field]

    # Note: Hydro Routing `SubId` is a float for some reason and Python float != GeoServer double. Cast them to int.
    if isinstance(fid, float):
        fid = int(fid)
        df[basin_field] = df[basin_field].astype(int)
        df[downstream_field] = df[downstream_field].astype(int)

    # Locate the downstream feature
    ds = df.set_index(basin_field).loc[fid]
    if basin_family is not None:
        # Do a first selection on the main basin ID of the downstream feature.
        sub = df[df[basin_family] == ds[basin_family]]
    else:
        sub = None

    # Find upstream basins
    up = [fid]
    for b in up:
        tmp = upstream_ids(sub if sub is not None else df, b)
        if len(tmp):
            up.extend(tmp)

    return (
        sub[sub[basin_field].isin(up)]
        if sub is not None
        else df[df[basin_field].isin(up)]
    )


def find_geometry_from_coord(
    lon: float, lat: float, df: geopandas.GeoDataFrame
) -> geopandas.GeoDataFrame:
    """
    Return the geometry containing the given coordinates.

    Parameters
    ----------
    lon : float
        Longitude.
    lat : float
        Latitude.
    df : GeoDataFrame
        Data.

    Returns
    -------
    geopandas.GeoDataFrame
        Record whose geometry contains the point.
    """

    p = Point(lon, lat)

    c = df.contains(p)
    n = c.sum()
    if n == 0:
        raise ValueError(f"Point {p} not in any geometry.")
    elif n > 1:
        raise ValueError(f"Point {p} found in multiple geometries.")

    return df[c]
