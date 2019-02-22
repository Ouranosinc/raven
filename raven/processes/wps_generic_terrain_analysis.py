import json
import logging
import tempfile
import rasterio
import geopandas as gpd
import rasterio.mask
import rasterio.warp
import shapely.ops as ops

from pywps import LiteralInput, ComplexInput, ComplexOutput
from pywps import Process, FORMATS
from pywps.app.Common import Metadata
from rasterio.crs import CRS
from raven.utils import archive_sniffer, crs_sniffer, single_file_check
from raven.utils import gdal_aspect_analysis, gdal_slope_analysis
from raven.utils import geom_transform
from shapely.geometry import shape


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
                         'Coordinate Reference System for terrain analysis (EPSG code; Default:32198).'
                         ' The CRS chosen should be projected and appropriate for the region of interest.',
                         data_type='integer',
                         default=32198,
                         min_occurs=1, max_occurs=1)
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
        shape_url = request.inputs['shape'][0].file
        destination_crs = request.inputs['projected_crs'][0].data

        rasters = ['.tiff', '.tif']
        raster_file = single_file_check(archive_sniffer(raster_url, working_dir=self.workdir, extensions=rasters))
        ras_crs = crs_sniffer(raster_file)

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

        if ras_crs == projection.to_proj4():
            reproject_raster = False
            msg = 'CRS for raster matches projected CRS {}. Raster will not be warped.'.format(projection.to_epsg())
            LOGGER.info(msg)
        else:
            reproject_raster = True
            msg = 'Warping raster to destination CRS {}'.format(projection.to_epsg())
            LOGGER.info(msg)

        if shape_url:
            clipped_raster = tempfile.NamedTemporaryFile(suffix='.tiff', delete=False)

            reproject_shape = False
            vectors = ['.gml', '.shp', '.geojson', '.json']  # '.gpkg' requires more handling
            vector_file = archive_sniffer(shape_url, working_dir=self.workdir, extensions=vectors)
            vec_crs = crs_sniffer(vector_file)
            if vec_crs != projection:
                reproject_shape = True

            if reproject_shape:
                # Reproject with GEOM functions
                pass

            gdf = gpd.GeoDataFrame.from_file(vector_file, crs=vec_crs)
            geom = shape(ops.unary_union(gdf['geometry']))
            transformed = False
            try:
                if reproject_shape:
                    x_min, y_min, x_max, y_max = geom.bounds
                    if y_max > 60 or y_min < 60:
                        # see https://gis.stackexchange.com/a/145138/65343 for discussion on this
                        msg = 'Shape boundaries exceed Mercator region.' \
                              'Verify choice of projected CRS is appropriate for aspect analysis.'
                        LOGGER.warning(msg)
                        raise UserWarning(msg)
                    transformed = geom_transform(geom, source_crs=vec_crs, target_crs=destination_crs)

                with rasterio.open(raster_file, 'r') as src:
                    if reproject_raster:
                        affine, width, height = rasterio.warp.calculate_default_transform(
                            ras_crs, projection, src.width, src.height, *src.bounds
                        )
                        mask_meta = src.meta.copy
                        mask_meta.update(
                            {
                                "driver": "GTiff",
                                "height": height,
                                "width": width,
                                "transform": affine
                            }
                        )

                        with rasterio.open(clipped_raster.name, 'w', **mask_meta) as dst:
                            rasterio.warp.reproject(
                                source=rasterio.band(src, 1),
                                destination=rasterio.band(dst, 1),
                                src_transform=src.transform,
                                src_crs=src.crs,
                                dst_transform=affine,
                                dst_crs=projection,
                                resampling=rasterio.warp.Resampling.nearest
                            )

                    else:
                        mask_image, mask_affine = rasterio.mask.mask(src, transformed or geom, crop=True)
                        mask_meta = src.meta.copy
                        mask_meta.update(
                            {
                                "driver": "GTiff",
                                "height": mask_image.shape[1],
                                "width": mask_image.shape[2],
                                "transform": mask_affine
                             }
                        )

                        with rasterio.open(clipped_raster.name, 'w', **mask_meta) as dst:
                            dst.write(mask_image)

            except Exception as e:
                msg = '{}: Failed to mask DEM {} with shape {}'.format(raster_file, vector_file, e)
                LOGGER.error(msg)

        slope_raster = tempfile.NamedTemporaryFile(suffix='.tiff', delete=False)
        gdal_slope_analysis(clipped_raster or raster_file, slope_raster, units='degree', band=1)
        aspect_raster = tempfile.NamedTemporaryFile(suffix='.tiff', delete=False)
        gdal_aspect_analysis(clipped_raster or raster_file, aspect_raster, zeroForFlat=False, band=1)

        response.outputs['slope'].data = json.dumps(slope_raster)
        response.outputs['aspect'].data = json.dumps(aspect_raster)

        return response
