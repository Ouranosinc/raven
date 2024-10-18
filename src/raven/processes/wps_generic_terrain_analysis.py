import json
import logging
import tempfile

import shapely.geometry as sgeo
import shapely.ops as ops
from pyproj.crs import CRS
from pywps import FORMATS, ComplexOutput, LiteralInput, Process

from raven.utilities.analysis import dem_prop
from raven.utilities.checks import boundary_check, single_file_check
from raven.utilities.geo import (
    generic_raster_clip,
    generic_raster_warp,
    generic_vector_reproject,
)
from raven.utilities.io import archive_sniffer, crs_sniffer
from raven.utils import gather_dem_tile

from . import wpsio as wio

LOGGER = logging.getLogger("PYWPS")


class TerrainAnalysisProcess(Process):
    """Given a file containing vector data and a DEM, analyze terrain characteristics."""

    def __init__(self):
        inputs = [
            wio.dem_raster,
            wio.shape,
            LiteralInput(
                "projected_crs",
                "Coordinate Reference System for terrain analysis (Default: EPSG:6622, 'NAD83(CSRS) /"
                " Quebec Lambert'. The CRS chosen should be projected and appropriate for the region"
                " of interest.",
                data_type="integer",
                default=6622,
                min_occurs=1,
                max_occurs=1,
            ),
            wio.select_all_touching,
        ]

        outputs = [
            ComplexOutput(
                "properties",
                "Feature schemas",
                abstract="DEM properties (mean elevation, slope, and aspect) for each geometry.",
                supported_formats=[FORMATS.JSON],
            ),
            ComplexOutput(
                "dem",
                "Subsetted digital elevation model",
                abstract="DEM GeoTIFF image",
                as_reference=True,
                supported_formats=[FORMATS.GEOTIFF, FORMATS.META4],
            ),
        ]

        super().__init__(
            self._handler,
            identifier="terrain-analysis",
            title="Terrain Analysis",
            version="1.0",
            abstract="Return shape area in square metres based on line boundaries of a polygonal vector file.",
            metadata=[],
            inputs=inputs,
            outputs=outputs,
            status_supported=True,
            store_supported=True,
        )

    def _handler(self, request, response):
        # Process inputs
        # ---------------
        shape_url = request.inputs["shape"][0].file
        destination_crs = request.inputs["projected_crs"][0].data
        touches = request.inputs["select_all_touching"][0].data

        # Checks for valid CRS and that CRS is projected
        # -----------------------------------------------
        projection = CRS.from_user_input(destination_crs)
        if not projection.is_projected:
            msg = f"Destination CRS {projection.to_epsg()} is not projected. Terrain analysis values will not be valid."
            LOGGER.error(ValueError(msg))
            raise ValueError(msg)

        # Collect and process the shape
        # -----------------------------
        vectors = [".gml", ".shp", ".gpkg", ".geojson", ".json"]
        vector_file = single_file_check(
            archive_sniffer(shape_url, working_dir=self.workdir, extensions=vectors)
        )
        vec_crs = crs_sniffer(vector_file)

        # Check that boundaries within 60N and 60S
        boundary_check(vector_file)

        if "raster" in request.inputs:
            raster_url = request.inputs["raster"][0].file
            rasters = [".tiff", ".tif"]
            raster_file = single_file_check(
                archive_sniffer(
                    raster_url, working_dir=self.workdir, extensions=rasters
                )
            )

        else:
            # Assuming that the shape coordinate are in WGS84
            raster_file = gather_dem_tile(vector_file, self.workdir)

        ras_crs = crs_sniffer(raster_file)

        # Reproject raster
        # ----------------
        if ras_crs != projection.to_epsg():
            msg = f"CRS for {raster_file} is not {projection}. Reprojecting raster..."
            LOGGER.warning(msg)
            warped_fn = tempfile.NamedTemporaryFile(
                prefix="warped_", suffix=".tiff", delete=False, dir=self.workdir
            ).name
            generic_raster_warp(raster_file, warped_fn, projection)

        else:
            warped_fn = raster_file

        # Perform the terrain analysis
        # ----------------------------
        rpj = tempfile.NamedTemporaryFile(
            prefix="reproj_", suffix=".json", delete=False, dir=self.workdir
        ).name

        generic_vector_reproject(
            vector_file, rpj, source_crs=vec_crs, target_crs=projection
        )
        with open(rpj) as src:
            geo = json.load(src)

        features = [sgeo.shape(feat["geometry"]) for feat in geo["features"]]
        union = ops.unary_union(features)

        clipped_fn = tempfile.NamedTemporaryFile(
            prefix="clipped_", suffix=".tiff", delete=False, dir=self.workdir
        ).name
        # Ensure that values for regions outside of clip are kept
        generic_raster_clip(
            raster=warped_fn,
            output=clipped_fn,
            geometry=union,
            touches=touches,
            fill_with_nodata=True,
            padded=True,
        )

        # Compute DEM properties for each feature.
        properties = []
        for i in range(len(features)):
            properties.append(
                dem_prop(clipped_fn, geom=features[i], directory=self.workdir)
            )
        properties.append(dem_prop(clipped_fn, directory=self.workdir))

        response.outputs["properties"].data = json.dumps(properties)
        response.outputs["dem"].file = clipped_fn

        return response
