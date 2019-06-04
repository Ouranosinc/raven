import json
import logging
import math
import os
import re
import tarfile
import tempfile
import zipfile
from functools import partial
from re import search

import fiona
import numpy as np
import pyproj
import rasterio
import rasterio.mask
import rasterio.vrt
import rasterio.warp
import shapely.geometry as sgeo
from osgeo.gdal import DEMProcessing
from rasterio.crs import CRS
from shapely.ops import transform

LOGGER = logging.getLogger("RAVEN")

# See: https://kokoalberti.com/articles/geotiff-compression-optimization-guide/
GDAL_TIFF_COMPRESSION_OPTION = 'compress=lzw'  # or 'compress=deflate' or 'compress=zstd' or 'compress=lerc' or others
RASTERIO_TIFF_COMPRESSION = 'lzw'

WGS84 = '+init=epsg:4326'
WGS84_PROJ4 = '+proj=longlat +datum=WGS84 +no_defs'
LAEA = '+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
WORLDMOLL = '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
ALBERS_NAM = '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs'


def address_append(address):
    """
    Formats a URL/URI to be more easily read with libraries such as "rasterstats"

    Parameters
    ----------
    address: str or path
      URL/URI to a potential zip or tar file

    Returns
    -------
    address: str
      URL/URI prefixed for archive type
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
    """Extract archives (tar/zip) to a working directory.

    Parameters
    ----------
    resources: list
      list of archive files (if netCDF files are in list, they are passed and returned as well in the return).
    output_dir: str or path
      string or Path to a working location (default: temporary folder).

    Returns
    -------
    list
      List of original or of extracted files
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
                    msg = '7z file extraction is not supported at this time'
                    LOGGER.warning(msg)
                    raise UserWarning(msg)
                else:
                    LOGGER.debug('File extension "{}" unknown'.format(file))
            except Exception as e:
                LOGGER.error('Failed to extract sub archive {}: {}'.format(arch, e))
        else:
            LOGGER.warning('No archives found. Continuing...')
            return resources

    return files


def archive_sniffer(archives, working_dir, extensions):
    """Return a list of locally unarchived files that match the desired extensions.

    Parameters
    ----------
    archives : str or path or list
      archive location or list of archive locations
    working_dir : str or path
      string or Path to a working location
    extensions : list
      list of accepted extensions

    Returns
    -------
    potential_files : list
      List of files with matching accepted extensions
    """
    potential_files = []

    decompressed_files = generic_extract_archive(archives, output_dir=working_dir)
    for file in decompressed_files:
        if any(ext in os.path.splitext(file) for ext in extensions):
            potential_files.append(file)
    return potential_files


def crs_sniffer(*args):
    """Return the list of CRS found in files.

    Parameters
    ----------
    args : sequence
      Paths to the files to examine.

    Returns
    -------
    crs_list : list or str
      Returns either a list of CRSes or a single CRS definition, depending on the number of instances found.
    """
    crs_list = []
    vectors = ('.gml', '.shp', '.geojson', '.gpkg', '.json')
    rasters = ('.tif', '.tiff')

    for file in args:
        found_crs = False
        try:
            if str(file).lower().endswith(vectors):
                if str(file).lower().endswith('.gpkg'):
                    if len(fiona.listlayers(file)) > 1:
                        msg = 'Multilayer GeoPackages are currently unsupported'
                        LOGGER.warning(msg)
                        raise NotImplementedError(msg)
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
        if crs_list[0] == '':
            msg = 'No CRS definitions found in {}. Using {}'.format(args, WGS84)
            LOGGER.warning(msg)
            return WGS84
        return crs_list[0]
    return crs_list


def raster_datatype_sniffer(file):
    """Return the type of the raster stored in the file.

    Parameters
    ----------
    file : str
      Path to file.

    Returns
    -------
    dtype: str
      rasterio datatype of array values
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
    """Return longitude and latitude from a string."""
    try:
        lon, lat = tuple(map(float, re.findall(r'[-+]?[0-9]*\.?[0-9]+', string)))
        return lon, lat
    except Exception as e:
        msg = '{}: Failed to parse longitude, latitude coordinates {}'.format(e, string)
        raise ValueError(msg)


def single_file_check(file_list):
    """Return the first element of a file list. Raise an error if the list is empty or contains more than one element.
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


def boundary_check(*args, max_y=60, min_y=-60):
    """Verify that boundaries do not exceed specific latitudes for geographic coordinate data. Raise a warning if so.

    Parameters
    ----------
    *args : str or path
      listing of strings or paths to files
    max_y : int or float
      Maximum value allowed for latitude. Default: 60.
    min_y : int or float
      Minimum value allowed for latitude. Default: -60.
    """
    vectors = ('.gml', '.shp', '.geojson', '.gpkg', '.json')
    rasters = ('.tif', '.tiff')
    for file in args:
        src = None
        try:
            if str(file).lower().endswith(vectors):
                src = fiona.open(file, 'r')
            elif str(file).lower().endswith(rasters):
                src = rasterio.open(file, 'r')
            else:
                FileNotFoundError('Invalid filename suffix')

            geographic = CRS(src.crs).is_geographic
            if geographic and (src.bounds > max_y or src.bounds < min_y):
                msg = 'Vector {} contains geometries in high latitudes.' \
                      ' Verify choice of projected CRS is appropriate for analysis.'.format(file)
                LOGGER.warning(msg)
                UserWarning(msg)
            src.close()

        except Exception as e:
            msg = '{}: Unable to read boundaries from {}'.format(e, args)
            LOGGER.exception(msg)
    return


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
def geom_transform(geom, source_crs=WGS84, target_crs=None):
    """Change the projection of a geometry.

    Assuming a geometry's coordinates are in a `source_crs`, compute the new coordinates under the `target_crs`.

    Parameters
    ----------
    geom : shapely.geometry
      Source geometry.
    source_crs : str or CRS
      Projection identifier (proj4) for the source geometry, e.g. '+proj=longlat +datum=WGS84 +no_defs'.
    target_crs : str or CRS
      Projection identifier (proj4) for the target geometry.

    Returns
    -------
    shapely.geometry
      Reprojected geometry.
    """
    try:
        reprojected = transform(
            partial(
                pyproj.transform,
                pyproj.Proj(source_crs),
                pyproj.Proj(target_crs)),
            geom)
        return reprojected
    except Exception as e:
        msg = '{}: Failed to reproject geometry'.format(e)
        LOGGER.error(msg)
        raise Exception(msg)


def geom_prop(geom):
    """Return a dictionary of geometry properties.

    Parameters
    ----------
    geom : shapely.geometry
      Geometry to analyze.

    Returns
    -------
    dict
      Dictionary storing polygon area, centroid location, perimeter and gravelius shape index.

    Notes
    -----
    Some of the properties should be computed using an equal-area projection.
    """

    geom = sgeo.shape(geom)
    lon, lat = geom.centroid.x, geom.centroid.y
    if (lon > 180) or (lon < -180) or (lat > 90) or (lat < -90):
        LOGGER.warning('Shape centroid is not in decimal degrees.')
    area = geom.area
    length = geom.length
    gravelius = length / 2 / math.sqrt(math.pi * area)
    parameters = {'area': area,
                  'centroid': (lon, lat),
                  'perimeter': length,
                  'gravelius': gravelius,
                  }
    return parameters


def dem_prop(dem, geom=None, directory=None):
    """Return raster properties for each geometry.

    This

    Parameters
    ----------
    dem : str
      DEM raster in reprojected coordinates.
    geom : shapely.geometry
      Geometry over which aggregate properties will be computed. If None compute properties over entire raster.
    directory : str or path
      Folder to save the GDAL terrain anlaysis outputs

    Returns
    -------
    dict
      Dictionary storing mean elevation [m], slope [deg] and aspect [deg].
    """

    fns = dict()
    fns['dem'] = tempfile.NamedTemporaryFile(prefix='dem', suffix='.tiff', dir=directory, delete=False).name \
        if geom is not None else dem
    for key in ['slope', 'aspect']:
        fns[key] = tempfile.NamedTemporaryFile(prefix=key, suffix='.tiff', dir=directory, delete=False).name

    # Clip to relevant area or read original raster
    if geom is None:
        with rasterio.open(dem) as f:
            elevation = f.read(1, masked=True)
    else:
        generic_raster_clip(raster=dem, output=fns['dem'], geometry=geom)
        with rasterio.open(fns['dem']) as f:
            elevation = f.read(1, masked=True)

    # Compute slope
    slope = gdal_slope_analysis(fns['dem'], output=fns['slope'])

    # Compute aspect
    aspect = gdal_aspect_analysis(fns['dem'], output=fns['aspect'])
    aspect_mean = circular_mean_aspect(aspect)

    return {'elevation': elevation.mean(), 'slope': slope.mean(), 'aspect': aspect_mean}


# It's a bit weird to have to pass the output file name as an argument, since you return an in-memory array.
# Can you keep everything in memory ?

# Response: The DEMPROCESSING function is basically a thin GDAL wrapper for a os.subprocess call. There is
# unfortunately no way to call the gdal slope/aspect calculation from Python directly so it demands an output filename.
# Since it technically is writing the information to the file, this function could use a generic named temporary file
# to perform the analysis and return the array all in memory I suppose. I've added this as an option here.

# @multipolygon_check
def gdal_slope_analysis(dem, output=None, units='degree'):
    """Return the slope of the terrain from the DEM.

    The slope is the magnitude of the gradient of the elevation.

    Parameters
    ----------
    dem : str
      Path to file storing DEM.
    output : str
      Path to output file.
    units : str
      Slope units. Default: 'degree'.

    Returns
    -------
    ndarray
      Slope array.

    Notes
    -----
    Ensure that the DEM is in a *projected coordinate*, not a geographic coordinate system, so that the
    horizontal scale is the same as the vertical scale (m).

    """
    if output is None:
        output = tempfile.NamedTemporaryFile().name
    DEMProcessing(output, dem, 'slope', slopeFormat=units,
                  format='GTiff', band=1, creationOptions=[GDAL_TIFF_COMPRESSION_OPTION, ])
    with rasterio.open(output) as src:
        return np.ma.masked_values(src.read(1), value=-9999)


def gdal_aspect_analysis(dem, output=None, flat_values_are_zero=False):
    """Return the aspect of the terrain from the DEM.

    The aspect is the compass direction of the steepest slope (0: North, 90: East, 180: South, 270: West).

    Parameters
    ----------
    dem : str
      Path to file storing DEM.
    output : str
      Path to output file.
    flat_values_are_zero: bool
      Designate flat values with value zero. Default: -9999.

    Returns
    -------
    ndarray
      Aspect array.

    Notes
    -----
    Ensure that the DEM is in a *projected coordinate*, not a geographic coordinate system, so that the
    horizontal scale is the same as the vertical scale (m).
    """
    if output is None:
        output = tempfile.NamedTemporaryFile().name
    DEMProcessing(destName=output, srcDS=dem, processing='aspect', zeroForFlat=flat_values_are_zero,
                  format='GTiff', band=1, creationOptions=[GDAL_TIFF_COMPRESSION_OPTION, ])
    with rasterio.open(output) as src:
        return np.ma.masked_values(src.read(1), value=-9999)


def circular_mean_aspect(angles):
    """Return the mean angular aspect based on circular arithmetic approach

    Parameters
    ----------
    angles: ndarray
      Array of aspect angles

    Returns
    -------
    float
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


def generic_raster_clip(raster, output, geometry, touches=False, fill_with_nodata=True, padded=True,
                        raster_compression=RASTERIO_TIFF_COMPRESSION):
    """
    Crop a raster file to a given geometry.

    Parameters
    ----------
    raster : str
      Path to input raster.
    output : str
      Path to output raster.
    geometry : shapely.geometry
      Geometry defining the region to crop.
    touches : bool
      Whether or not to include cells that intersect the geometry. Default: True.
    fill_with_nodata: bool
      Whether or not to keep pixel values for regions outside of shape or set as nodata. Default: True.
    padded: bool
      Whether or not to add a half-pixel buffer to shape before masking
    raster_compression : str
      Level of data compression. Default: 'lzw'.

    """
    if not (type(geometry) in (list, tuple)):
        geometry = [geometry]

    with rasterio.open(raster, 'r', ) as src:
        mask_image, mask_affine = rasterio.mask.mask(src, geometry, crop=True, pad=padded,
                                                     all_touched=touches, filled=fill_with_nodata)
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
        with rasterio.open(output, 'w', **mask_meta) as dst:
            dst.write(mask_image)
    return


def generic_raster_warp(raster, output, target_crs, raster_compression=RASTERIO_TIFF_COMPRESSION):
    """
    Reproject a raster file.

    Parameters
    ----------
    raster : str
      Path to input raster.
    output : str
      Path to output raster.
    target_crs : str or dict
      Target projection identifier.
    raster_compression: str
      Level of data compression. Default: 'lzw'.
    """
    with rasterio.open(raster, 'r') as src:
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

            with rasterio.open(output, 'w', **metadata) as dst:
                dst.write(data)
    return


def generic_vector_reproject(vector, projected, source_crs=WGS84, target_crs=None):
    """Reproject all features and layers within a vector file and return a GeoJSON

    Parameters
    ----------
    vector : str or path
      Path to a file containing a valid vector layer.
    projected: str or path
      Path to a file to be written.
    source_crs : str or CRS
      Projection identifier (proj4) for the source geometry, Default: '+proj=longlat +datum=WGS84 +no_defs'.
    target_crs : str or CRS
      Projection identifier (proj4) for the target geometry.

    Returns
    -------
    None
    """

    if target_crs is None:
        msg = 'No target CRS is defined.'
        raise ValueError(msg)

    output = {"type": "FeatureCollection", 'features': []}

    for i, layer_name in enumerate(fiona.listlayers(vector)):
        with fiona.open(vector, 'r', layer=i) as src:
            with open(projected, 'w') as sink:
                for feature in src:
                    # Perform vector reprojection using Shapely on each feature
                    try:
                        geom = sgeo.shape(feature['geometry'])
                        transformed = geom_transform(geom, source_crs, target_crs)
                        feature['geometry'] = sgeo.mapping(transformed)
                        output['features'].append(feature)
                    except Exception as e:
                        msg = '{}: Unable to reproject feature {}'.format(e, feature)
                        LOGGER.exception(msg)

                sink.write("{}".format(json.dumps(output)))
    return
