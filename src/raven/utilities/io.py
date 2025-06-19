"""Tools for reading and writing geospatial data formats."""

import logging
import os
import tarfile
import tempfile
import warnings
import zipfile
from collections.abc import Sequence
from pathlib import Path
from re import search
from typing import Optional, Union

import fiona
import rasterio
from pyproj import CRS
from shapely.geometry import shape

LOGGER = logging.getLogger("RavenPy")
WGS84 = 4326


# Function addressing exploit CVE-2007-4559
def is_within_directory(
    directory: Union[str, os.PathLike], target: Union[str, os.PathLike]
) -> bool:
    abs_directory = os.path.abspath(directory)
    abs_target = os.path.abspath(target)

    prefix = os.path.commonprefix([abs_directory, abs_target])

    return prefix == abs_directory


# Function addressing exploit CVE-2007-4559
def safe_extract(
    tar: tarfile.TarFile, path: str = ".", members=None, *, numeric_owner=False
) -> None:
    for member in tar.getmembers():
        member_path = os.path.join(path, member.name)
        if not is_within_directory(path, member_path):
            raise Exception("Attempted Path Traversal in Tar File")

    tar.extractall(path, members, numeric_owner=numeric_owner)


def address_append(address: Union[str, Path]) -> str:
    """
    Format a URL/URI to be more easily read with libraries such as "rasterstats".

    Parameters
    ----------
    address : Union[str, Path]
        URL/URI to a potential zip or tar file.

    Returns
    -------
    str
        URL/URI prefixed for the archive type.
    """
    zipped = search(r"(\.zip)", str(address))
    tarred = search(r"(\.tar)", str(address))

    try:
        if zipped:
            return f"zip://{address}"
        elif tarred:
            return f"tar://{address}"
        else:
            LOGGER.info("No prefixes needed for address.")
            return str(address)
    except Exception:
        LOGGER.error("Failed to prefix or parse URL %s." % address)
        raise


def generic_extract_archive(
    resources: Union[str, Path, list[Union[bytes, str, Path]]],
    output_dir: Optional[Union[str, Path]] = None,
) -> list[str]:
    """
    Extract archives (tar/zip) to a working directory.

    Parameters
    ----------
    resources : str or Path or list of bytes or str or Path
        List of archive files (if netCDF files are in a list, they are passed and returned as well in the return).
    output_dir : str or Path, optional
        String or Path to a working location (default: temporary folder).

    Returns
    -------
    list
        A list of original or of extracted files.
    """

    archive_types = [".tar", ".zip", ".7z"]
    output_dir = output_dir or tempfile.gettempdir()

    if not isinstance(resources, list):
        resources = [resources]

    files = list()

    for arch in resources:
        if any(ext in str(arch).lower() for ext in archive_types):
            try:
                LOGGER.debug("archive=%s", arch)
                file = Path(arch).name

                if file.endswith(".nc"):
                    files.append(Path(output_dir.join(arch)))
                elif file.endswith(".tar"):
                    with tarfile.open(arch, mode="r") as tar:
                        safe_extract(
                            tar, path=output_dir
                        )  # Function addressing exploit CVE-2007-4559
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
                    msg = "7z file extraction is not supported at this time."
                    LOGGER.warning(msg)
                    warnings.warn(msg, UserWarning)
                else:
                    LOGGER.debug('File extension "%s" unknown' % file)
            except Exception as e:
                LOGGER.error(f"Failed to extract sub archive {{{arch}}}: {{{e}}}")
        else:
            LOGGER.warning("No archives found. Continuing...")
            return resources

    return files


def archive_sniffer(
    archives: Union[str, Path, list[Union[str, Path]]],
    working_dir: Optional[Union[str, Path]] = None,
    extensions: Optional[Sequence[str]] = None,
) -> list[Union[str, Path]]:
    """
    Return a list of locally unarchived files that match the desired extensions.

    Parameters
    ----------
    archives : str or Path or list of str or Path
        Archive location or list of archive locations.
    working_dir : str or Path, optional
        String or Path to a working location.
    extensions : Sequence of str, optional
        List of accepted extensions.

    Returns
    -------
    list of str or Path
        List of files with matching accepted extensions.
    """
    potential_files = list()

    if not extensions:
        extensions = [".gml", ".shp", ".geojson", ".gpkg", ".json"]

    decompressed_files = generic_extract_archive(archives, output_dir=working_dir)
    for file in decompressed_files:
        if any(ext in Path(file).suffix for ext in extensions):
            potential_files.append(file)
    return potential_files


def crs_sniffer(
    *args: Union[str, Path, Sequence[Union[str, Path]]],
) -> Union[list[Union[str, int]], str, int]:
    """Return the list of CRS found in files.

    Parameters
    ----------
    *args : Union[str, Path, Sequence[Union[str, Path]]]
        Path(s) to the file(s) to examine.

    Returns
    -------
    Union[List[str], str]
        Returns either a list of CRSes or a single CRS definition, depending on the number of instances found.
    """
    crs_list = list()
    vectors = (".gml", ".shp", ".geojson", ".gpkg", ".json")
    rasters = (".tif", ".tiff")
    all_files = vectors + rasters

    for file in args:
        found_crs = False
        suffix = Path(file).suffix.lower()
        try:
            if suffix == ".zip":
                file = archive_sniffer(file, extensions=all_files)[0]
                suffix = Path(file).suffix.lower()

            if suffix in vectors:
                if suffix == ".gpkg":
                    if len(fiona.listlayers(file)) > 1:
                        raise NotImplementedError
                with fiona.open(file, "r") as src:
                    found_crs = CRS.from_wkt(src.crs_wkt).to_epsg()
            elif suffix in rasters:
                with rasterio.open(file, "r") as src:
                    found_crs = CRS.from_user_input(src.crs).to_epsg()
            else:
                raise FileNotFoundError("Invalid filename suffix")
        except FileNotFoundError as e:
            msg = f"{e}: Unable to open file {args}"
            LOGGER.warning(msg)
            raise Exception(msg)
        except NotImplementedError as e:
            msg = f"{e}: Multilayer GeoPackages are currently unsupported"
            LOGGER.error(msg)
            raise Exception(msg)
        except RuntimeError:
            pass

        crs_list.append(found_crs)

    if crs_list is None:
        msg = f"No CRS definitions found in {args}."
        raise FileNotFoundError(msg)

    if len(crs_list) == 1:
        if not crs_list[0]:
            msg = f"No CRS definitions found in {args}. Assuming {WGS84}."
            LOGGER.warning(msg)
            warnings.warn(msg, UserWarning)
            return WGS84
        return crs_list[0]
    return crs_list


def raster_datatype_sniffer(file: Union[str, Path]) -> str:
    """
    Return the type of the raster stored in the file.

    Parameters
    ----------
    file : Union[str, Path]
        Path to file.

    Returns
    -------
    str
        The rasterio datatype of array values.
    """
    try:
        with rasterio.open(file, "r") as src:
            dtype = src.dtypes[0]
        return dtype
    except rasterio.errors.RasterioError:
        msg = f"Unable to read data type from {file}."
        LOGGER.exception(msg)
        raise ValueError(msg)


def get_bbox(
    vector: Union[str, Path], all_features: bool = True
) -> tuple[float, float, float, float]:
    """
    Return bounding box of all features or the first feature in file.

    Parameters
    ----------
    vector : str or Path
        A path to file storing vector features.
    all_features : bool
        Return the bounding box for all features. Default: True.

    Returns
    -------
    float, float, float, float
        Geographic coordinates of the bounding box (lon0, lat0, lon1, lat1).
    """

    if not all_features:
        with fiona.open(vector, "r") as src:
            for feature in src:
                geom = shape(feature["geometry"])
                return geom.bounds

    with fiona.open(vector, "r") as src:
        return src.bounds
