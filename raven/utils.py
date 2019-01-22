from re import search
from functools import partial

import tempfile
import logging
import math
import os
import tarfile
import zipfile

import shapely.geometry as sgeo
import shapely.ops as ops
from fiona.crs import from_epsg
import pyproj

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


def extract_archive(resources, output_dir=None):
    """
    extracts archives (tar/zip)
    :param resources: list of archive files (if netCDF files are in list,
                     they are passed and returned as well in the return).
    :param output_dir: define a directory to store the results (default: tempory folder).
    :return list: [list of extracted files]
    """
    output_dir = output_dir or tempfile.gettempdir()

    if not isinstance(resources, list):
        resources = list([resources])
    files = []

    for arch in resources:
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
            else:
                LOGGER.warning('file extension "{}" unknown'.format(file))
        except Exception as e:
            LOGGER.error('failed to extract sub archive {}: {}'.format(arch, e))
    return files


def archive_sniffer(file, working_dir, archive_types, extensions):
    if any(ext in str(file).lower() for ext in archive_types):
        extracted = extract_archive(file, working_dir)
        for potential_file in extracted:
            if any(str(potential_file).lower().endswith(ext) for ext in extensions):
                return potential_file


def multipolygon_check(geom):
    """Perform a check to verify a geometry is a MultiPolygon

    Parameters
    ----------
    geom : shapely.geometry
    """
    if isinstance(sgeo.shape(geom), sgeo.multipolygon.MultiPolygon):
        LOGGER.warning("Shape is a Multipolygon.")
    return


def geom_prop(geom):
    """Return a dictionary of properties for the given geometry.

    Parameters
    ----------
    geom : shapely.geometry
    """
    multipolygon_check(geom)

    geom = sgeo.shape(geom)
    out = {'centroid': (geom.centroid.x, geom.centroid.y)}
    return out


def geom_transform(geom, source_crs, target_crs):
    """Return a projected geometry based on source and target CRS

    Parameters
    ----------
    geom : shapely.geometry
    source_crs : EPSG code (Default: 4326)
    target_crs : EPSG code

    Returns
    -------
    shapely.geometry
      A projected shapely.geometry
    """
    if not source_crs:
        source_crs = 4326

    projected = ops.transform(
        partial(
            pyproj.transform,
            pyproj.Proj(from_epsg(source_crs)),
            pyproj.Proj(from_epsg(target_crs))),
        geom)

    return projected


def equal_area_geom_prop(geom):
    """Return a dictionary of properties for the given equal area geometry.

    Parameters
    ----------
    geom : shapely.geometry

    Returns
    -------
    dict
      Dictionary storing polygon area, perimeter and gravelius shape index.
    """
    multipolygon_check(geom)

    geom = sgeo.shape(geom)
    area = geom.area
    length = geom.length
    gravelius = length / 2 / math.sqrt(math.pi * area)
    parameters = {'area': area,
                  'perimeter': length,
                  'gravelius': gravelius,
                  }
    return parameters
