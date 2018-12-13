from re import search
import tempfile
import logging
import math
import os
from zipfile import ZipFile
import tarfile

LOGGER = logging.getLogger("RAVEN")


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
                with ZipFile(arch, mode='r') as zf:
                    zf.extractall(path=output_dir)
                    files.extend([os.path.join(output_dir, f) for f in zf.namelist()])
            else:
                LOGGER.warning('file extension "{}" unknown'.format(file))
        except Exception as e:
            LOGGER.error('failed to extract sub archive {}: {}'.format(arch, e))
    return files


def geom_prop(geom):
    """Return a dictionary of properties for the given geometry.

    Parameters
    ----------
    geom : shapely.geometry
    """
    # TODO: Check that it's not a MULTIPOLYGON
    out = {'centroid': (geom.centroid.x, geom.centroid.y)}
    return out


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
    # TODO: Check that it's not a MULTIPOLYGON
    area = geom.area
    length = geom.length
    gravelius = length / 2 / math.sqrt(math.pi * area)
    out = {'area': area,
           'perimeter': length,
           'gravelius': gravelius,
           }
    return out
