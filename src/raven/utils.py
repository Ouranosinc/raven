import logging
import os
import re
import shutil
from pathlib import Path
from random import choice
from string import ascii_letters
from tempfile import NamedTemporaryFile
from typing import Union

import numpy as np
import rasterio
import rasterio.mask
import rasterio.vrt
import rasterio.warp
from affine import Affine
from pyproj.crs import CRS

from raven.utilities.geoserver import get_raster_wcs
from raven.utilities.io import get_bbox

LOGGER = logging.getLogger("RAVEN")

RASTERIO_TIFF_COMPRESSION = "lzw"

TRUE_CATEGORIES = {
    0: "Ocean",
    1: "Temperate or sub-polar needleleaf forest",
    2: "Sub-polar taiga needleleaf forest",
    3: "Tropical or sub-tropical broadleaf evergreen forest",
    4: "Tropical or sub-tropical broadleaf deciduous forest",
    5: "Temperate or sub-polar broadleaf deciduous forest",
    6: "Mixed forest",
    7: "Tropical or sub-tropical shrubland",
    8: "Temperate or sub-polar shrubland",
    9: "Tropical or sub-tropical grassland",
    10: "Temperate or sub-polar grassland",
    11: "Sub-polar or polar shrubland-lichen-moss",
    12: "Sub-polar or polar grassland-lichen-moss",
    13: "Sub-polar or polar barren-lichen-moss",
    14: "Wetland",
    15: "Cropland",
    16: "Barren lands",
    17: "Urban",
    18: "Water",
    19: "Snow and Ice",
}

simplified = {
    "Ocean": [0],
    "Forest": [1, 2, 3, 4, 5, 6],
    "Shrubs": [7, 8, 11],
    "Grass": [9, 10, 12, 13, 16],
    "Wetland": [14],
    "Crops": [15],
    "Urban": [17],
    "Water": [18],
    "SnowIce": [19],
}
SIMPLE_CATEGORIES = {i: cat for (cat, ids) in simplified.items() for i in ids}

SUMMARY_ZONAL_STATS = ["count", "nodata", "nan"]
NALCMS_PROJ4 = (
    "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs=True"
)
EARTH_ENV_DEM = "public:EarthEnv_DEM90_NorthAmerica"


def gather_dem_tile(
    vector_file: Union[str, os.PathLike],
    work_dir: Union[str, os.PathLike],
    geographic: bool = True,
    raster: str = EARTH_ENV_DEM,
) -> Path:
    """Return a raster coverage for a given vector geometry.

    Parameters
    ----------
    vector_file : str or os.PathLike
        Shape whose bounds will be used to collect a raster coverage.
    work_dir : str or os.PathLike
        Folder where the file will be written.
    geographic : bool
        Use geographic units (degree-decimal) or projected units (metres/feet).
    raster : str
        Layer name on GeoServer.

    Returns
    -------
    Path
        Path to raster file.
    """
    bbox = get_bbox(vector_file)
    raster_layer = raster
    raster_bytes = get_raster_wcs(bbox, geographic=geographic, layer=raster_layer)
    raster_file = NamedTemporaryFile(
        prefix="wcs_", suffix=".tiff", delete=False, dir=work_dir
    ).name
    with open(raster_file, "wb") as f:
        f.write(raster_bytes)
    return Path(raster_file)


def parse_lonlat(lonlat: Union[str, tuple[str, str]]) -> tuple[float, float]:
    """Return longitude and latitude from a string.

    Parameters
    ----------
    lonlat : str or Tuple[str, str]
      A tuple or a str of lon and lat coordinates.

    Returns
    -------
    (float, float)
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
        msg = f"Failed to parse longitude, latitude coordinates {lonlat}"
        raise Exception(msg) from e


def zonalstats_raster_file(
    stats: dict,
    working_dir: str = None,
    raster_compression: str = RASTERIO_TIFF_COMPRESSION,
    data_type: str = None,
    crs: str = None,
    zip_archive: bool = False,
) -> Union[Path, list[Path]]:
    """
    Extract the zonalstats grid(s) to a zipped GeoTIFF file and ensure that it is projected with specified CRS.

    Parameters
    ----------
    stats : dict
      The dictionary produced by the rasterstats `zonalstats` function.
    working_dir : str
      The working directory.
    raster_compression : str
      The type of compression used on the raster file (default: 'lzw').
    data_type : str
      The data encoding of the raster used to write the grid (e.g. 'int16').
    crs : str
      The coordinate reference system.
    zip_archive : bool
      Whether to return the files as a zipped archive (default: False).

    Returns
    -------
    Path or List[Path]
    """
    out_dir = Path(working_dir).joinpath("output")
    out_dir.mkdir(exist_ok=True)
    crs = CRS(crs)

    for i in range(len(stats)):
        fn = f"subset_{i + 1}.tiff"
        raster_subset = Path(out_dir).joinpath(fn)

        try:
            raster_location = stats[i]
            raster = raster_location["mini_raster_array"]
            grid_properties = raster_location["mini_raster_affine"][0:6]
            nodata = raster_location["mini_raster_nodata"]

            aff = Affine(*grid_properties)

            LOGGER.info(f"Writing raster data to {raster_subset}")

            masked_array = np.ma.masked_values(raster, nodata)
            if masked_array.mask.all():
                msg = f"Subset {i} is empty, continuing..."
                LOGGER.warning(msg)

            normal_array = np.asarray(masked_array, dtype=data_type)

            # Write to GeoTIFF
            with rasterio.open(
                raster_subset,
                "w",
                driver="GTiff",
                count=1,
                compress=raster_compression,
                height=raster.shape[0],
                width=raster.shape[1],
                dtype=data_type,
                transform=aff,
                crs=crs,
                nodata=nodata,
            ) as f:
                f.write(normal_array, 1)

        except Exception as e:
            msg = f"Failed to write raster outputs: {e}"
            LOGGER.error(msg)
            raise Exception(msg)

    # `shutil.make_archive` could potentially cause problems with multi-thread? Worth investigating later.
    if zip_archive:
        foldername = f"subset_{''.join(choice(ascii_letters) for _ in range(10))}"
        out_fn = Path(working_dir).joinpath(foldername)
        shutil.make_archive(
            base_name=out_fn.as_posix(), format="zip", root_dir=out_dir, logger=LOGGER
        )
        return Path(f"{out_fn}.zip")
    else:
        return [f for f in out_dir.glob("*")]
