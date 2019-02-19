from re import search
from functools import partial

import tempfile
import logging
import math
import os
import re
import tarfile
import zipfile

import shapely.geometry as sgeo
import shapely.ops as ops
import fiona as fio
from fiona.crs import from_epsg
from pyproj import Proj, transform
import rasterio as rio

LOGGER = logging.getLogger("RAVEN")

LAEA = '+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
WORLDMOLL = '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
ALBERS_NAM = '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs'


def address_append(address):
    """
    Formats a URL/URI to be more easily read with libraries such as "rasterstats"
    :param address: URL/URI to a potential zip or tar file
    :return: URL/URI prefixed for archive type
    """
    zipped = search(r'(\.zip)', address)
    tarred = search(r'(\.tar)', address)

    try:
        if zipped:
            return 'zip://{}'.format(address)
        elif tarred:
            return 'tar://{}'.format(address)
        else:
            LOGGER.info('No changes made to address.')
            return address
    except Exception as e:
        msg = 'Failed to prefix or parse URL {}: {}'.format(address, e)
        LOGGER.error(msg)


def generic_extract_archive(resources, output_dir=None):
    """
    Extracts archives (tar/zip) to a working directory
    :param resources: list of archive files (if netCDF files are in list,
                     they are passed and returned as well in the return).
    :param output_dir: string or Path to a working location (default: temporary folder).
    :return list: [list of extracted files]
    """
    archive_types = ['.tar', '.zip', '.7z']
    output_dir = output_dir or tempfile.gettempdir()

    if not isinstance(resources, list):
        resources = list([resources])

    files = []

    for arch in resources:
        if any(ext in str(arch).lower() for ext in archive_types):
            try:
                LOGGER.debug("archive=%s", arch)
                file = os.path.basename(arch)

                if file.endswith('.nc'):
                    files.append(os.path.join(output_dir, arch))
                elif file.endswith('.tar'):
                    with tarfile.open(arch, mode='r') as tar:
                        tar.extractall(path=output_dir)
                        files.extend([os.path.join(output_dir, f) for f in tar.getnames()])
                elif file.endswith('.zip'):
                    with zipfile.ZipFile(arch, mode='r') as zf:
                        zf.extractall(path=output_dir)
                        files.extend([os.path.join(output_dir, f) for f in zf.namelist()])
                elif file.endswith('.7z'):
                    LOGGER.warning('7z file extraction is not supported at this time')
                else:
                    LOGGER.debug('File extension "{}" unknown'.format(file))
            except Exception as e:
                LOGGER.error('Failed to extract sub archive {}: {}'.format(arch, e))
        else:
            LOGGER.warning('No archives found. Continuing...')
            return resources

    return files


def archive_sniffer(archives, working_dir, extensions):
    """
    Return a list of locally unarchived files that match the desired extensions
    :param archives : archive location or list of archive locations
    :param working_dir: string or Path to a working location
    :param extensions: [list of accepted extensions]
    :return:
    """
    potential_files = []

    decompressed_files = generic_extract_archive(archives, output_dir=working_dir)
    for file in decompressed_files:
        if any(ext in os.path.splitext(file) for ext in extensions):
            return file
    return potential_files


def crs_sniffer(file):
    found_crs = False

    vectors = ('.gml', '.shp', '.geojson', '.json')
    rasters = ('.tif', '.tiff')

    try:
        if file.lower().endswith(vectors):
            with fio.open(file, 'r') as src:
                found_crs = src.crs
                src.close()
        elif file.lower().endswith(rasters):
            with rio.open(file, 'r') as src:
                found_crs = src.crs
                src.close()
        else:
            raise FileNotFoundError('Invalid file suffix')
    except Exception as e:
        msg = '{}: Unable to read crs from {}'.format(e, file)
        LOGGER.warning(msg)
        return
    finally:
        if found_crs:  # Empty strings are treated as False
            return found_crs


def dtype_sniffer(file):
    try:
        with rio.open(file, 'r') as src:
            dtype = src.dtypes
        return dtype
    except Exception as e:
        msg = '{}: Unable to read data type from {}'.format(e, file)
        LOGGER.warning(msg)
        return


def parse_lonlat(string):
    try:
        lon, lat = tuple(map(float, re.findall(r'[-+]?[0-9]*\.?[0-9]+', string)))
        return lon, lat
    except Exception as e:
        msg = 'Failed to parse geo-coordinates {}: {}'.format(string, e)
        raise Exception(msg)


def multipolygon_check(f):
    pass
    # def wrapper(*args, **kwargs):
    #     """Perform a check to verify a geometry is a MultiPolygon
    #
    #        :params *args: shapely.geometry
    #        """
    #     if isinstance(type(*args), sgeo.multipolygon.MultiPolygon):
    #         LOGGER.warning("Shape is a Multipolygon.")
    #     print('ok!')
    # return wrapper


# @multipolygon_check
def geom_prop(geom):
    """
    Return a dictionary of properties for the given geometry.
    :param geom : shapely.geometry
    :return dict : centroid with lon, lat fields
    """
    shape = sgeo.shape(geom)
    out = {'centroid': (shape.centroid.x, shape.centroid.y)}
    if (out['centroid'][0] > 180) or (out['centroid'][0] < -180)\
            or (out['centroid'][1] > 90) or (out['centroid'][1] < -90):
        LOGGER.warning('Shape centroid is not in decimal degrees.')
    return out


# @multipolygon_check
def geom_transform(geom, source_crs=4326, target_crs=None):
    """
    Return a projected geometry based on source and target CRS
    :param geom : shapely.geometry
    :param source_crs : EPSG code (Default: 4326)
    :param target_crs : EPSG code
    :return shapely.geometry:
    """
    geom = sgeo.shape(geom)
    projected = ops.transform(
        partial(
            transform,
            Proj(from_epsg(source_crs)),
            Proj(from_epsg(target_crs))),
        geom)
    return projected


# @multipolygon_check
def equal_area_geom_prop(geom):
    """
    Return a dictionary of properties for the given equal area geometry.
    :param geom : shapely.geometry
    :return dict : Dictionary storing polygon area, perimeter and gravelius shape index.
    """
    geom = sgeo.shape(geom)
    area = geom.area
    length = geom.length
    gravelius = length / 2 / math.sqrt(math.pi * area)
    parameters = {'area': area,
                  'perimeter': length,
                  'gravelius': gravelius,
                  }
    return parameters
