import json
import logging
import os
import tempfile
import shutil

import numpy as np
import rasterio as rio

from affine import Affine
from pywps import LiteralInput, ComplexInput
from pywps import ComplexOutput
from pywps import Process, FORMATS
from pywps.app.Common import Metadata
from rasterstats import zonal_stats

from raven.utils import archive_sniffer, crs_sniffer, single_file_check, raster_datatype_sniffer

LOGGER = logging.getLogger("PYWPS")


class RasterSubsetProcess(Process):
    """Given files containing vector data and raster data, perform zonal statistics of the overlapping regions."""

    def __init__(self):
        inputs = [
            ComplexInput('shape', 'Vector Shape',
                         abstract='An ESRI Shapefile, GML, GeoJSON, or any other file in a standard vector format.'
                                  ' The ESRI Shapefile must be zipped and contain the .shp, .shx, and .dbf.'
                                  ' The shape and raster should have a matching CRS.',
                         min_occurs=1, max_occurs=1,
                         supported_formats=[FORMATS.GEOJSON, FORMATS.GML, FORMATS.JSON, FORMATS.SHP]),
            ComplexInput('raster', 'Gridded raster data set',
                         abstract='The DEM to be subset. Defaults to the USGS HydroSHEDS DEM.',
                         # TODO: Include details (resolution, version).
                         metadata=[Metadata('HydroSheds Database', 'http://hydrosheds.org'),
                                   Metadata(
                                       'Lehner, B., Verdin, K., Jarvis, A. (2008): New global hydrography derived from'
                                       ' spaceborne elevation data. Eos, Transactions, AGU, 89(10): 93-94.',
                                       'https://doi.org/10.1029/2008EO100001')],
                         min_occurs=1, max_occurs=1, supported_formats=[FORMATS.GEOTIFF]),
            LiteralInput('band', 'Raster band', data_type='integer', default=1,
                         abstract='Band of raster examined to perform zonal statistics. Default: 1',
                         min_occurs=1, max_occurs=1),
            LiteralInput('select_all_touching', 'Additionally select boundary pixels that are touched by shape',
                         data_type='boolean', default='false'),
        ]

        outputs = [
            ComplexOutput('raster', 'DEM subset of `shape` region in GeoTIFF format.',
                          abstract='Elevation statistics: min, max, mean, median, sum, nodata',
                          supported_formats=[FORMATS.ZIP, ]),
        ]

        super(RasterSubsetProcess, self).__init__(
            self._handler,
            identifier="raster-subset",
            title="Raster Subset",
            version="1.0",
            abstract="Return a masked raster based on boundaries of a vector file.",
            metadata=[],
            inputs=inputs,
            outputs=outputs,
            status_supported=True,
            store_supported=True)

    def _handler(self, request, response):

        shape_url = request.inputs['shape'][0].file
        raster_url = request.inputs['raster'][0].file
        band = request.inputs['band'][0].data
        touches = request.inputs['select_all_touching'][0].data

        vectors = ['.gml', '.shp', '.geojson', '.json']  # '.gpkg' requires more handling
        rasters = ['.tiff', '.tif']
        vector_file = single_file_check(archive_sniffer(shape_url, working_dir=self.workdir, extensions=vectors))
        raster_file = single_file_check(archive_sniffer(raster_url, working_dir=self.workdir, extensions=rasters))

        vec_crs, ras_crs = crs_sniffer(vector_file), crs_sniffer(raster_file)

        if ras_crs != vec_crs:
            msg = 'CRS for files {} and {} are not the same.'.format(vector_file, raster_file)
            LOGGER.warning(msg)

        data_type = raster_datatype_sniffer(raster_file)
        raster_compression = 'lzw'

        out_dir = os.path.join(self.workdir, 'output')
        os.makedirs(out_dir)

        try:
            stats = zonal_stats(
                vector_file, raster_file, band=band, all_touched=touches, raster_out=True)

            for i in range(len(stats)):

                file = 'subset_{}.tiff'.format(i + 1)
                raster_subset = os.path.join(out_dir, file)

                try:
                    raster_location = stats[i]
                    raster = raster_location['mini_raster_array']
                    grid_properties = raster_location['mini_raster_affine'][0:6]
                    nodata = raster_location['mini_raster_nodata']

                    aff = Affine(*grid_properties)

                    LOGGER.info('Writing raster data to {}'.format(raster_subset))

                    masked_array = np.ma.masked_values(raster, nodata)
                    if masked_array.mask.all():
                        msg = 'Subset {} is empty, continuing...'.format(i)
                        LOGGER.warning(msg)

                    normal_array = np.asarray(masked_array, dtype=data_type)

                    # Write to GeoTIFF
                    with rio.open(raster_subset, 'w', driver='GTiff', count=1, compress=raster_compression,
                                  height=raster.shape[0], width=raster.shape[1], dtype=data_type, transform=aff,
                                  crs=vec_crs or ras_crs, nodata=nodata) as f:
                        f.write(normal_array, 1)

                except Exception as e:
                    msg = 'Failed to write raster outputs: {}'.format(e)
                    LOGGER.error(msg)
                    raise Exception(msg)

            # `shutil.make_archive` could potentially cause problems with multi-thread? Worth investigating later.
            out_fn = os.path.join(self.workdir, self.identifier)
            shutil.make_archive(base_name=out_fn, format='zip', root_dir=out_dir, logger=LOGGER)

            response.outputs['raster'].file = out_fn + '.zip'

        except Exception as e:
            msg = 'Failed to perform raster subset using {} and {}: {}'.format(shape_url, raster_url, e)
            LOGGER.error(msg)
            raise Exception(msg)

        return response
