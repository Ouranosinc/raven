import json
import logging
import tempfile

from pywps import LiteralInput, Process
from rasterstats import zonal_stats

from raven.utilities.checks import single_file_check
from raven.utilities.geo import generic_vector_reproject
from raven.utilities.io import archive_sniffer, crs_sniffer
from raven.utils import gather_dem_tile

from . import wpsio as wio

LOGGER = logging.getLogger("PYWPS")
SUMMARY_ZONAL_STATS = ["count", "min", "max", "mean", "median", "sum", "nodata"]


class ZonalStatisticsProcess(Process):
    """Given files containing vector data and raster data, perform zonal statistics of the overlapping regions"""

    def __init__(self):
        inputs = [
            wio.shape,
            wio.dem_raster,
            wio.raster_band,
            LiteralInput(
                "categorical",
                "Return distinct pixel categories",
                data_type="boolean",
                default="false",
                min_occurs=1,
                max_occurs=1,
            ),
            wio.select_all_touching,
        ]

        outputs = [
            wio.statistics,
        ]

        super().__init__(
            self._handler,
            identifier="zonal-stats",
            title="Raster Zonal Statistics",
            version="1.0",
            abstract="Return zonal statistics based on the boundaries of a vector file.",
            metadata=[],
            inputs=inputs,
            outputs=outputs,
            status_supported=True,
            store_supported=True,
        )

    def _handler(self, request, response):
        shape_url = request.inputs["shape"][0].file
        band = request.inputs["band"][0].data
        categorical = request.inputs["categorical"][0].data
        touches = request.inputs["select_all_touching"][0].data

        vectors = [".gml", ".shp", ".gpkg", ".geojson", ".json"]
        vector_file = single_file_check(
            archive_sniffer(shape_url, working_dir=self.workdir, extensions=vectors)
        )
        rasters = [".tiff", ".tif"]

        if "raster" in request.inputs:
            raster_url = request.inputs["raster"][0].file
            raster_file = single_file_check(
                archive_sniffer(
                    raster_url, working_dir=self.workdir, extensions=rasters
                )
            )
        else:
            raster_url = None
            # Assuming that the shape coordinate are in WGS84
            raster_file = gather_dem_tile(vector_file, self.workdir, geographic=True)

        vec_crs, ras_crs = crs_sniffer(vector_file), crs_sniffer(raster_file)

        if ras_crs != vec_crs:
            msg = f"CRS for files {vector_file} and {raster_file} are not the same. Reprojecting vector..."
            LOGGER.warning(msg)

            # Reproject full vector to preserve feature attributes
            projected = tempfile.NamedTemporaryFile(
                prefix="reprojected_", suffix=".json", delete=False, dir=self.workdir
            ).name
            generic_vector_reproject(
                vector_file, projected, source_crs=vec_crs, target_crs=ras_crs
            )
            vector_file = projected

        summary_stats = SUMMARY_ZONAL_STATS

        try:
            stats = zonal_stats(
                vector_file,
                raster_file,
                stats=summary_stats,
                band=band,
                categorical=categorical,
                all_touched=touches,
                geojson_out=True,
                raster_out=False,
            )

            feature_collect = {"type": "FeatureCollection", "features": stats}
            response.outputs["statistics"].data = json.dumps(feature_collect)

        except Exception as e:
            msg = f"Failed to perform raster subset using {shape_url}{f' and {raster_url} ' if raster_url else ''}: {e}"
            LOGGER.error(msg)
            raise Exception(msg) from e

        return response
