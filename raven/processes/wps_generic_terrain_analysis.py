import json
import logging
import os
import tempfile

# import numpy as np
# import rasterio
import geopandas as gpd
# import rasterio.mask
# import rasterio.warp
import shapely.ops as ops
import shapely.geometry as sgeo

from pywps import LiteralInput, ComplexInput, ComplexOutput
from pywps import Process, FORMATS
from pywps.app.Common import Metadata
from rasterio.crs import CRS
from raven.utils import archive_sniffer, crs_sniffer, single_file_check
from raven.utils import gdal_aspect_analysis, gdal_slope_analysis
from raven.utils import generic_raster_warp, generic_raster_clip, dem_prop
# from shapely.geometry import shape


LOGGER = logging.getLogger("PYWPS")


class TerrainAnalysisProcess(Process):
    """Given a file containing vector data and a DEM, analyze terrain characteristics."""

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
                         'Coordinate Reference System for terrain analysis (Default: EPSG:32198,'
                         ' "NAD83 / Quebec Lambert")).'
                         ' The CRS chosen should be projected and appropriate for the region of interest.',
                         data_type='integer',
                         default=32198,
                         min_occurs=1, max_occurs=1),
            LiteralInput('select_all_touching', 'Perform calculation on boundary pixels',
                         data_type='boolean', default='false',
                         min_occurs=1, max_occurs=1),
        ]

        outputs = [
            ComplexOutput('properties', 'Feature schemas',
                          abstract='DEM properties (mean elevation, slope and aspect) for each geometry.',
                          supported_formats=[FORMATS.JSON],
                          ),
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

        # Process inputs
        # ---------------
        raster_url = request.inputs['raster'][0].file

        try:
            shape_url = request.inputs['shape'][0].file
        except KeyError:
            shape_url = False

        destination_crs = request.inputs['projected_crs'][0].data
        # touches = request.inputs['select_all_touching'][0].data

        # Checks for valid CRS and that CRS is projected
        projection = CRS.from_user_input(destination_crs)
        if not projection.is_projected:
            msg = 'Destination CRS {} is not projected.' \
                  'Terrain analysis values will not be valid.'.format(projection.to_epsg())
            raise ValueError(msg)

        # Get vector features
        # -------------------
        if shape_url:
            vectors = ['.gml', '.shp', '.geojson', '.json']  # '.gpkg' requires more handling
            vector_file = single_file_check(archive_sniffer(shape_url, working_dir=self.workdir, extensions=vectors))
            vec_crs = crs_sniffer(vector_file)[0]

            # Load shape data with GeoPandas
            gdf = gpd.GeoDataFrame.from_file(vector_file, crs=vec_crs)

            # TODO: This check seems specific to a given projection, no ?
            x_min, y_min, x_max, y_max = sgeo.shape(ops.unary_union(gdf['geometry'])).bounds
            if y_max > 60 or y_min < 60:
                # see https://gis.stackexchange.com/a/145138/65343 for discussion on this
                msg = 'Shape boundaries exceed Mercator region.' \
                      'Verify choice of projected CRS is appropriate for aspect analysis.'
                LOGGER.warning(msg)
                UserWarning(msg)

            # Reproject geometries.
            if vec_crs != projection.to_proj4():
                geoms = gdf.to_crs(epsg=destination_crs)
            else:
                geoms = gdf

            union = sgeo.shape(ops.unary_union(geoms['geometry']))

        else:
            geoms = [None, ]

        # Raster handling
        # ---------------

        # Verify that the CRS to project to matches shape or raster CRS
        rasters = ['.tiff', '.tif']
        raster_file = single_file_check(archive_sniffer(raster_url, working_dir=self.workdir, extensions=rasters))
        ras_crs = crs_sniffer(raster_file)

        # Clip raster to unioned geometry
        if shape_url:
            clipped_fn = tempfile.NamedTemporaryFile(prefix='clipped_', suffix='.tiff', delete=False,
                                                     dir=self.workdir).name
            generic_raster_clip(raster_file, clipped_fn, union)

        else:
            clipped_fn = raster_file

        # Reproject raster
        if ras_crs != projection.to_proj4():
            # processed_raster = os.path.join(self.workdir, 'warped_{}.tiff'.format(hash(os.times())))
            warped_fn = tempfile.NamedTemporaryFile(prefix='warped_', suffix='.tiff', delete=False,
                                                    dir=self.workdir).name
            generic_raster_warp(clipped_fn, warped_fn, projection)
        else:
            warped_fn = clipped_fn

        # Compute DEM properties for each feature.
        properties = []
        for geom in geoms:
            properties.append(dem_prop(warped_fn, geom))

        response.outputs['properties'].data = json.dumps(properties)

        return response
