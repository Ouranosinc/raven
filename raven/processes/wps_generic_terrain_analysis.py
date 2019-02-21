# import fiona
import json
import logging

from rasterio.crs import CRS
from pywps import LiteralInput, ComplexInput, ComplexOutput
from pywps import Process, FORMATS
from pywps.app.Common import Metadata
from raven.utils import archive_sniffer, crs_sniffer, single_file_check
# from raven.utils import geom_transform, geom_centroid, equal_area_geom_prop
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
                         'Coordinate Reference System for terrain analysis (EPSG code; Default:32198).'
                         ' The CRS chosen should be projected and appropriate for the region of interest.',
                         data_type='integer',
                         default=32198,
                         min_occurs=1, max_occurs=1)
        ]

        outputs = [
            ComplexOutput('terrain_analysis', 'Terrain analysis characteristics of the DEM',
                          abstract='Terrain analysis characteristics of the DEM (Slope, Aspect, and Curvature)',
                          supported_formats=[FORMATS.JSON]),
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
            msg = 'CRS for raster matches projected CRS ({}). Raster will not be warped.'.format(projection.to_epsg())
            LOGGER.info(msg)
        else:
            msg = 'Warping raster to destination CRS ({})'.format(projection.to_epsg())
            LOGGER.info(msg)

        if shape_url:
            reproject_shape = False
            vectors = ['.gml', '.shp', '.geojson', '.json']  # '.gpkg' requires more handling
            vector_file = archive_sniffer(shape_url, working_dir=self.workdir, extensions=vectors)
            vec_crs = crs_sniffer(vector_file)
            if vec_crs != projection:
                reproject_shape = True

            if reproject_shape:
                # Reproject with GEOM functions
                pass

            properties = [].append(vec_crs)

            # try:
            #     with fiona.open(vector_file, 'r', crs=vec_crs) as src:
            #         for feature in src:
            #             geom = shape(feature['geometry'])
            #
            #             transformed = geom_transform(geom, source_crs=vec_crs, target_crs=destination_crs)
            #             prop = {'id': feature['id']}
            #             prop.update(feature['properties'])
            #             prop.update(geom_centroid(geom))
            #             prop.update(equal_area_geom_prop(transformed))
            #             properties.append(prop)
            #             break
            #
            # except Exception as e:
            #     msg = 'Failed to extract shape from url {}: {}'.format(vector_file, e)
            #     LOGGER.error(msg)

            response.outputs['terrain_analysis'].data = json.dumps(properties)

        return response
