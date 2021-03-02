import logging
import shutil
from pathlib import Path
from random import choice
from string import ascii_letters
from typing import List, Union

import numpy as np
import rasterio
import rasterio.mask
import rasterio.vrt
import rasterio.warp
from affine import Affine
from pyproj.crs import CRS

LOGGER = logging.getLogger("RAVEN")

RASTERIO_TIFF_COMPRESSION = "lzw"


def zonalstats_raster_file(
    stats: dict,
    working_dir: str = None,
    raster_compression: str = RASTERIO_TIFF_COMPRESSION,
    data_type: str = None,
    crs: str = None,
    zip_archive: bool = False,
) -> Union[str, List[Path]]:
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
    zip_archive: bool
      Return the files as a zipped archive (default: False).

    Returns
    -------
    Union[str, List[Path]]
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
            base_name=out_fn, format="zip", root_dir=out_dir, logger=LOGGER
        )
        return f"{out_fn}.zip"
    else:
        return [f for f in out_dir.glob("*")]
