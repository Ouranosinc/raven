import json
import logging
import tempfile

import shapely.geometry as sgeo
import shapely.ops as ops
from pywps import LiteralInput, ComplexInput, ComplexOutput
from pywps import Process, FORMATS
from pywps.app.Common import Metadata
from rasterio.crs import CRS

from raven.utilities import gis
from raven.utils import archive_sniffer, crs_sniffer, single_file_check, boundary_check
from raven.utils import generic_raster_warp, generic_raster_clip, dem_prop, generic_vector_reproject

LOGGER = logging.getLogger("PYWPS")


class TerrainAnalysisProcess(Process):
    """Given a file containing vector data and a DEM, analyze terrain characteristics."""

    def __init__(self):
        inputs = [
            ComplexInput('raster', 'Digital elevation model (DEM)',
                         abstract='The DEM to be queried. Defaults to the EarthEnv-DEM90 product.',
                         metadata=[Metadata('EarthEnv-DEM90', 'https://www.earthenv.org/DEM'),
                                   Metadata(
                                       'Robinson, Natalie, James Regetz, and Robert P. Guralnick (2014). '
                                       'EarthEnv-DEM90: A Nearly-Global, Void-Free, Multi-Scale Smoothed, 90m Digital '
                                       'Elevation Model from Fused ASTER and SRTM Data. ISPRS Journal of '
                                       'Photogrammetry and Remote Sensing 87: 57–67.',
                                       'https://doi.org/10.1016/j.isprsjprs.2013.11.002')],
                         min_occurs=0, max_occurs=1, supported_formats=[FORMATS.GEOTIFF]),
            ComplexInput('shape', 'Vector Shape',
                         abstract='An ESRI Shapefile, GML, JSON, GeoJSON, or single layer GeoPackage.'
                                  ' The ESRI Shapefile must be zipped and contain the .shp, .shx, and .dbf.'
                                  ' The raster will be subsetted before analysis is performed.',
                         min_occurs=1, max_occurs=1,
                         supported_formats=[FORMATS.GEOJSON, FORMATS.GML, FORMATS.JSON, FORMATS.SHP]),
            LiteralInput('projected_crs',
                         'Coordinate Reference System for terrain analysis (Default: EPSG:6622, "NAD83(CSRS) /'
                         ' Quebec Lambert". The CRS chosen should be projected and appropriate for the region'
                         ' of interest.',
                         data_type='integer',
                         default=6622,
                         min_occurs=1, max_occurs=1),
            LiteralInput('select_all_touching', 'Perform calculation on boundary pixels',
                         data_type='boolean', default='false',
                         min_occurs=1, max_occurs=1),
        ]

        outputs = [
            ComplexOutput('properties', 'Feature schemas',
                          abstract='DEM properties (mean elevation, slope, and aspect) for each geometry.',
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
        shape_url = request.inputs['shape'][0].file
        destination_crs = request.inputs['projected_crs'][0].data
        touches = request.inputs['select_all_touching'][0].data

        # Checks for valid CRS and that CRS is projected
        # -----------------------------------------------
        projection = CRS.from_user_input(destination_crs)
        if not projection.is_projected:
            msg = 'Destination CRS {} is not projected.' \
                  ' Terrain analysis values will not be valid.'.format(projection.to_epsg())
            LOGGER.error(ValueError(msg))
            raise ValueError(msg)

        # Collect and process the shape
        # -----------------------------
        vectors = ['.gml', '.shp', '.gpkg', '.geojson', '.json']
        vector_file = single_file_check(archive_sniffer(shape_url, working_dir=self.workdir, extensions=vectors))
        vec_crs = crs_sniffer(vector_file)

        # Check that boundaries within 60N and 60S
        boundary_check(vector_file)

        if 'raster' in request.inputs:
            raster_url = request.inputs['raster'][0].file
            rasters = ['.tiff', '.tif']
            raster_file = single_file_check(archive_sniffer(raster_url, working_dir=self.workdir, extensions=rasters))

        else:
            # Assuming that the shape coordinate are in WGS84
            bbox = gis.get_bbox(vector_file)
            raster_url = 'public:EarthEnv_DEM90_NorthAmerica'
            raster_bytes = gis.get_raster_wcs(bbox, geographic=True, layer=raster_url)
            raster_file = tempfile.NamedTemporaryFile(prefix='wcs_', suffix='.tiff', delete=False,
                                                      dir=self.workdir).name
            with open(raster_file, 'wb') as f:
                f.write(raster_bytes)

        ras_crs = crs_sniffer(raster_file)

        # Reproject raster
        # ----------------
        if ras_crs != projection.to_proj4():
            msg = 'CRS for {} is not {}. Reprojecting raster...'.format(raster_file, projection)
            LOGGER.warning(msg)
            warped_fn = tempfile.NamedTemporaryFile(prefix='warped_', suffix='.tiff', delete=False,
                                                    dir=self.workdir).name
            generic_raster_warp(raster_file, warped_fn, projection.to_proj4())

        else:
            warped_fn = raster_file

        # Perform the terrain analysis
        # ----------------------------
        rpj = tempfile.NamedTemporaryFile(prefix='reproj_', suffix='.json', delete=False, dir=self.workdir).name
        generic_vector_reproject(vector_file, rpj, source_crs=vec_crs, target_crs=projection.to_proj4())
        with open(rpj) as src:
            geo = json.load(src)

        features = [sgeo.shape(feat['geometry']) for feat in geo['features']]
        union = ops.unary_union(features)

        clipped_fn = tempfile.NamedTemporaryFile(prefix='clipped_', suffix='.tiff', delete=False,
                                                 dir=self.workdir).name
        # Ensure that values for regions outside of clip are kept
        generic_raster_clip(raster=warped_fn, output=clipped_fn, geometry=features, touches=touches,
                            fill_with_nodata=True, padded=True)

        # Compute DEM properties for each feature.
        properties = []
        for i in range(len(features)):
            properties.append(dem_prop(clipped_fn, features[i], directory=self.workdir))
        properties.append(dem_prop(clipped_fn, directory=self.workdir))

        response.outputs['properties'].data = json.dumps(properties)

        return response
