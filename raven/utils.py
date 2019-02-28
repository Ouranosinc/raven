from re import search
from functools import partial

import tempfile
import logging
import math
import os
import re
import tarfile
import zipfile

import fiona
import rasterio
import rasterio.mask
import rasterio.warp
import shapely.geometry as sgeo
import shapely.ops as ops

from osgeo.gdal import DEMProcessing
from numpy import zeros_like
from pyproj import Proj, transform
from rasterio.crs import CRS

LOGGER = logging.getLogger("RAVEN")

# See: https://kokoalberti.com/articles/geotiff-compression-optimization-guide/
GDAL_TIFF_COMPRESSION_OPTION = 'compress=lzw'  # or 'compress=deflate' or 'compress=zstd' or 'compress=lerc' or others
RASTERIO_TIFF_COMPRESSION = 'lzw'

WGS84 = '+proj=longlat +datum=WGS84 +no_defs'
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
            potential_files.append(file)
    return potential_files


def crs_sniffer(*args):
    """

    :param args:
    :return:
    """
    crs_list = []
    vectors = ('.gml', '.shp', '.geojson', '.json')
    rasters = ('.tif', '.tiff')

    for file in args:
        found_crs = False
        try:
            if str(file).lower().endswith(vectors):
                with fiona.open(file, 'r') as src:
                    found_crs = CRS(src.crs).to_proj4()
                    src.close()
            elif str(file).lower().endswith(rasters):
                with rasterio.open(file, 'r') as src:
                    found_crs = CRS(src.crs).to_proj4()
                    src.close()
            else:
                FileNotFoundError('Invalid filename suffix')
        except Exception as e:
            msg = '{}: Unable to read crs from {}'.format(e, args)
            LOGGER.exception(msg)
        finally:
            crs_list.append(found_crs)

    if crs_list is None:
        msg = 'No CRS definitions found in {}.'.format(args)
        raise FileNotFoundError(msg)

    if len(crs_list) == 1:
        return crs_list[0]
    return crs_list


def raster_datatype_sniffer(file):
    """

    :param file:
    :return:
    """
    try:
        with rasterio.open(file, 'r') as src:
            dtype = src.dtypes[0]
        return dtype
    except Exception as e:
        msg = '{}: Unable to read data type from {}'.format(e, file)
        LOGGER.exception(msg)
        raise ValueError(msg)


def parse_lonlat(string):
    try:
        lon, lat = tuple(map(float, re.findall(r'[-+]?[0-9]*\.?[0-9]+', string)))
        return lon, lat
    except Exception as e:
        msg = '{}: Failed to parse LonLat-coordinates {}'.format(e, string)
        raise ValueError(msg)


def single_file_check(file_list):
    """

    :param file_list:
    :return:
    """
    try:
        if len(file_list) > 1:
            msg = 'Multi-file handling for file is not supported. Exiting.'
            raise NotImplementedError(msg)
        elif len(file_list) == 0:
            msg = 'No files found. Exiting.'
            raise FileNotFoundError(msg)
        return file_list[0]
    except Exception as e:
        msg = '{}: Unspecified error. Exiting,'
        LOGGER.error(msg.format(e))
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
def geom_centroid(geom):
    """
    Return a dictionary entry for centroid for the given geometry.
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
def geom_transform(geom, source_crs=WGS84, target_crs=None):
    """
    Return a projected geometry based on source and target CRS
    :param geom : shapely.geometry
    :param source_crs : EPSG code (Default: 4326)
    :param target_crs : EPSG code
    :return shapely.geometry:
    """
    geom = sgeo.shape(geom)

    try:
        projected = ops.transform(
            partial(
                transform,
                Proj(source_crs),
                Proj(target_crs)),
            geom)
    except Exception as e:
        msg = '{}: Failed to perform projection. Exiting,'.format(e)
        LOGGER.error(msg)
        raise Exception(msg)
    return projected


# @multipolygon_check
def geom_equal_area_prop(geom):
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


def gdal_slope_analysis(dem, output, units='degree'):
    """

    :param dem:
    :param output:
    :param units:
    :return:
    """
    DEMProcessing(output, dem, 'slope', slopeFormat=units,
                  format='GTiff', band=1, creationOptions=[GDAL_TIFF_COMPRESSION_OPTION, ])
    with rasterio.open(output) as src:
        slope = src.read(1)
    return slope


def gdal_aspect_analysis(dem, output, flat_values_are_zero=False):
    """
    :param dem:
    :param output:
    :param flat_values_are_zero:
    :return:
    """
    DEMProcessing(output, dem, 'aspect', zeroForFlat=flat_values_are_zero,
                  format='GTiff', band=1, creationOptions=[GDAL_TIFF_COMPRESSION_OPTION, ])
    with rasterio.open(output) as src:
        aspect = src.read(1)
    return aspect


def generic_raster_clip(raster_file, processed_raster, geometry, touches=True,
                        raster_compression=RASTERIO_TIFF_COMPRESSION):
    """

    :param raster_file:
    :param processed_raster:
    :param geometry:
    :param touches:
    :param raster_compression:
    :return:
    """
    with rasterio.open(raster_file, 'r') as src:

        mask_image, mask_affine = rasterio.mask.mask(src, geometry, crop=True,
                                                     all_touched=touches)
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
        with rasterio.open(processed_raster, 'w', **mask_meta) as dst:
            dst.write(mask_image)


def generic_raster_warp(raster_file, processed_raster, projection, raster_compression=RASTERIO_TIFF_COMPRESSION):
    """

    :param raster_file:
    :param processed_raster:
    :param projection:
    :param raster_compression:
    :return:
    """
    with rasterio.open(raster_file, 'r') as src:

        affine, width, height = rasterio.warp.calculate_default_transform(
            src.crs, projection, src.width, src.height, *src.bounds
        )
        metadata = src.meta.copy()
        metadata.update(
            {
                "driver": "GTiff",
                "height": height,
                "width": width,
                "transform": affine,
                "crs": projection,
                "compress": raster_compression,
            }
        )
        data = src.read()

        with rasterio.open(processed_raster, 'w', **metadata) as dst:
            for i, band in enumerate(data, 1):
                # Create an empty array
                dest = zeros_like(band)

                # Warp raster with crs/affine and fill nodata from source values
                rasterio.warp.reproject(
                    source=band,
                    destination=dest,
                    src_nodata=src.nodata,
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_nodata=src.nodata,
                    dst_transform=affine,
                    dst_crs=projection,
                    resampling=rasterio.warp.Resampling.nearest
                )

                dst.write(dest, indexes=i)
