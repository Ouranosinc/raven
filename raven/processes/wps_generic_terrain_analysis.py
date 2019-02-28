import json
import logging
import tempfile

# import numpy as np
# import rasterio
# import geopandas as gpd
# import rasterio.mask
# import rasterio.warp
# import shapely.ops as ops

from pywps import LiteralInput, ComplexInput, ComplexOutput
from pywps import Process, FORMATS
from pywps.app.Common import Metadata
from rasterio.crs import CRS
from raven.utils import archive_sniffer, crs_sniffer, single_file_check
from raven.utils import gdal_aspect_analysis, gdal_slope_analysis
from raven.utils import generic_raster_warp  # geom_transform, generic_raster_clip
# from shapely.geometry import shape


LOGGER = logging.getLogger("PYWPS")


class TerrainAnalysisProcess(Process):
    """Given a file containing vector data and a DEM, analyze terrain characteristics"""

    def __init__(self):
        inputs = [
            ComplexInput('raster', 'Digital elevation model (DEM)',
                         abstract='The DEM to be analyzed. Defaults to the USGS HydroSHEDS DEM.',
                         # TODO: Include details (resolution, version).
                         metadata=[Metadata('HydroSheds Database', 'http://hydrosheds.org'),
                                   Metadata(
                                       'Lehner, B., Verdin, K., Jarvis, A. (2008): New global hydrography derived from'
                                       ' spaceborne elevation data. Eos, Transactions, AGU, 89(10): 93-94.',
                                       'https://doi.org/10.1029/2008EO100001')],
                         min_occurs=1, max_occurs=1, supported_formats=[FORMATS.GEOTIFF]),
            ComplexInput('shape', 'Vector Shape',
                         abstract='An ESRI Shapefile, GML, GeoJSON, or any other file in a standard vector format.'
                                  ' The ESRI Shapefile must be zipped and contain the .shp, .shx, and .dbf.'
                                  ' If a shapefile is provided, the raster will be subsetted before analysis.',
                         min_occurs=0, max_occurs=1,
                         supported_formats=[FORMATS.GEOJSON, FORMATS.GML, FORMATS.JSON, FORMATS.SHP]),
            LiteralInput('projected_crs',
                         # TODO: Write the name of the EPSG CRS
                         'Coordinate Reference System for terrain analysis (EPSG code; Default:32198).'
                         ' The CRS chosen should be projected and appropriate for the region of interest.',
                         data_type='integer',
                         default=32198,
                         min_occurs=1, max_occurs=1),
            LiteralInput('select_all_touching', 'Perform calculation on boundary pixels',
                         data_type='boolean', default='false',
                         min_occurs=1, max_occurs=1),
        ]

        outputs = [
            ComplexOutput('slope', 'Slope values for the DEM',
                          abstract='Terrain analysis characteristics of the DEM (Slope, Aspect, and Curvature)',
                          supported_formats=[FORMATS.JSON]),
            ComplexOutput('aspect', 'Aspect values for the DEM',
                          abstract='Terrain analysis characteristics of the DEM (Slope, Aspect, and Curvature)',
                          supported_formats=[FORMATS.JSON])
        ]

        super(TerrainAnalysisProcess, self).__init__(
            self._handler,
            identifier="terrain-analysis",
            title="Terrain Analysis",
            version="1.0",
            abstract="Return shape area in square metres based on line boundaries of a polygonal vector file.",
            metadata=[],
            inputs=inputs,
            outputs=outputs,
            status_supported=True,
            store_supported=True)

    def _handler(self, request, response):

        raster_url = request.inputs['raster'][0].file
        try:
            shape_url = request.inputs['shape'][0].file
        except KeyError:
            shape_url = False
        destination_crs = request.inputs['projected_crs'][0].data
        # touches = request.inputs['select_all_touching'][0].data

        rasters = ['.tiff', '.tif']
        raster_file = single_file_check(archive_sniffer(raster_url, working_dir=self.workdir, extensions=rasters))
        ras_crs = crs_sniffer(raster_file)[0]
        raster_compression = 'lzw'

        processed_raster = False

        # Checks for valid CRS and that CRS is projected
        try:
            projection = CRS.from_user_input(destination_crs)

            if not projection.is_projected:
                msg = 'Destination CRS {} is not projected.' \
                      'Terrain analysis values may be erroneous.'.format(projection.to_epsg())
                LOGGER.warning(msg)
        except Exception as e:
            msg = '{}: Failed to parse CRS definition. Exiting.'.format(e)
            LOGGER.error(msg)
            raise Exception(msg)

        # Verify that the CRS to project to matches shape or raster CRS
        if ras_crs == projection.to_proj4():
            reproject_raster = False
            msg = 'CRS for raster matches projected CRS {}. Raster will not be warped.'.format(projection.to_epsg())
            LOGGER.info(msg)
        else:
            reproject_raster = True
            msg = 'Warping raster to destination CRS {}'.format(projection.to_epsg())
            LOGGER.info(msg)

        if shape_url or reproject_raster:
            processed_raster = tempfile.NamedTemporaryFile(suffix='.tiff', delete=False)

        if shape_url:
            msg = "This needs a bit more work but it's nearly finished"
            raise NotImplementedError(msg)

            # reproject_shape = False
            # vectors = ['.gml', '.shp', '.geojson', '.json']  # '.gpkg' requires more handling
            # vector_file = single_file_check(archive_sniffer(shape_url, working_dir=self.workdir, extensions=vectors))
            # vec_crs = crs_sniffer(vector_file)[0]
            #
            # if vec_crs != projection.to_proj4():
            #     reproject_shape = True
            #     msg = 'CRS for raster does not match projected CRS {}. Shape will be reprojected.'\
            #         .format(projection.to_epsg())
            #     LOGGER.info(msg)
            #
            # # Load shape data with GeoPandas
            # gdf = gpd.GeoDataFrame.from_file(vector_file, crs=vec_crs)
            # geometry = shape(ops.unary_union(gdf['geometry']))
            # transformed_shape = []
            #
            # try:
            #
            #     # Note that high latitudes require a particular CRS type to return sane values
            #     if reproject_shape:
            #         x_min, y_min, x_max, y_max = geometry.bounds
            #         if y_max > 60 or y_min < 60:
            #             # see https://gis.stackexchange.com/a/145138/65343 for discussion on this
            #             msg = 'Shape boundaries exceed Mercator region.' \
            #                   'Verify choice of projected CRS is appropriate for aspect analysis.'
            #             LOGGER.warning(msg)
            #             UserWarning(msg)
            #
            #         # Perform the shape reprojection
            #         transformed_shape.append(
            #             geom_transform(geometry, source_crs=vec_crs, target_crs=projection.to_proj4()))
            #
            #     with rasterio.open(raster_file, 'r') as src:
            #         if reproject_raster:
            #             # Warp the grid using the source and destination CRS definitions
            #             affine, width, height = rasterio.warp.calculate_default_transform(
            #                 ras_crs, projection, src.width, src.height, *src.bounds
            #             )
            #             mask_meta = src.meta.copy()
            #             mask_meta.update(
            #                 {
            #                     "driver": "GTiff",
            #                     "height": height,
            #                     "width": width,
            #                     "transform": affine
            #                 }
            #             )
            #
            #             # Warp and write new raster using the transformed grid
            #             with rasterio.open(processed_raster.name, 'w', **mask_meta) as dst:
            #                 rasterio.warp.reproject(
            #                     source=rasterio.band(src, band),
            #                     destination=rasterio.band(dst, band),
            #                     src_transform=src.transform,
            #                     src_crs=src.crs,
            #                     dst_transform=affine,
            #                     dst_crs=projection,
            #                     resampling=rasterio.warp.Resampling.nearest
            #                 )
            #                 mask_image, mask_affine = rasterio.mask.mask(
            #                 dst, transformed_shape or geometry, crop=True, all_touched=touches)
            #                 dst.write(mask_image)
            #
            #         # Only mask the raster image using the transformed or original shape
            #         if not reproject_raster:
            #             generic_raster_clip(raster_file, processed_raster.name, transformed_shape or geometry,
            #                                 raster_compression=raster_compression)
            #
            # except Exception as e:
            #     msg = '{}: Failed to clip and mask DEM {} with shape {}'.format(e, raster_file, vector_file)
            #     LOGGER.error(msg)
            #     raise Exception(msg)

        # If no shape to subset region, warp raster or provide link to raw raster file
        elif reproject_raster:
            try:
                generic_raster_warp(raster_file, processed_raster.name, projection,
                                    raster_compression=raster_compression)

            except Exception as e:
                msg = '{}: Failed to warp DEM {} to crs {}'.format(e, raster_file, projection.to_epsg())
                LOGGER.error(msg)
                raise Exception(msg)
        else:
            processed_raster.name = raster_file

        try:
            slope_raster = tempfile.NamedTemporaryFile(prefix='slope', suffix='.tiff', delete=False)
            gdal_slope_analysis(processed_raster.name, slope_raster.name, units='degree')
            aspect_raster = tempfile.NamedTemporaryFile(prefix='aspect', suffix='.tiff', delete=False)
            gdal_aspect_analysis(processed_raster.name, aspect_raster.name, flat_values_are_zero=False)
        except Exception as e:
            msg = '{}: Failed to calculate slope/aspect with  {} and crs {}'\
                .format(e, raster_file, projection.to_epsg())
            LOGGER.error(msg)
            raise Exception(msg)

        response.outputs['slope'].data = json.dumps(slope_raster.name)
        response.outputs['aspect'].data = json.dumps(aspect_raster.name)

        return response
