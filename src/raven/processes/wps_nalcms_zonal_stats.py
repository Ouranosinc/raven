import json
import logging
import tempfile
from collections import defaultdict

from pywps import Process
from rasterstats import zonal_stats

from raven.utilities import geoserver
from raven.utilities.checks import single_file_check
from raven.utilities.geo import generic_vector_reproject
from raven.utilities.io import archive_sniffer, crs_sniffer, get_bbox
from raven.utils import (
    NALCMS_PROJ4,
    SIMPLE_CATEGORIES,
    SUMMARY_ZONAL_STATS,
    TRUE_CATEGORIES,
)

from . import wpsio as wio

LOGGER = logging.getLogger("PYWPS")


class NALCMSZonalStatisticsProcess(Process):
    """Given files containing vector data and raster data, perform zonal statistics of the overlapping regions"""

    def __init__(self):
        inputs = [
            wio.shape,
            wio.land_use_raster,
            wio.simple_categories,
            wio.raster_band,
            wio.select_all_touching,
        ]

        outputs = [wio.features, wio.statistics]

        super().__init__(
            self._handler,
            identifier="nalcms-zonal-stats",
            title="NALCMS Land Use Zonal Statistics",
            version="1.0",
            abstract="Return zonal statistics and land-use cover for the CEC NALCMS based on the boundaries of a vector file.",
            metadata=[],
            inputs=inputs,
            outputs=outputs,
            status_supported=True,
            store_supported=True,
        )

    def _handler(self, request, response):
        shape_url = request.inputs["shape"][0].file
        simple_categories = request.inputs["simple_categories"][0].data
        band = request.inputs["band"][0].data
        touches = request.inputs["select_all_touching"][0].data

        vectors = [".gml", ".shp", ".gpkg", ".geojson", ".json"]
        vector_file = single_file_check(
            archive_sniffer(shape_url, working_dir=self.workdir, extensions=vectors)
        )
        vec_crs = crs_sniffer(vector_file)

        response.update_status("Accessed vector", status_percentage=5)

        if (
            "raster" in request.inputs
        ):  # For raster files using the UNFAO Land Cover Classification System (19 types)
            rasters = [".tiff", ".tif"]
            raster_url = request.inputs["raster"][0].file
            raster_file = single_file_check(
                archive_sniffer(
                    raster_url, working_dir=self.workdir, extensions=rasters
                )
            )
            ras_crs = crs_sniffer(raster_file)

            if vec_crs != ras_crs:
                msg = f"CRS for files {vector_file} and {raster_file} are not the same. Reprojecting..."
                LOGGER.warning(msg)

                # Reproject full vector to preserve feature attributes
                projected = tempfile.NamedTemporaryFile(
                    prefix="reprojected_",
                    suffix=".json",
                    delete=False,
                    dir=self.workdir,
                ).name
                generic_vector_reproject(
                    vector_file, projected, source_crs=vec_crs, target_crs=ras_crs
                )
            else:
                projected = vector_file

        else:  # using the NALCMS data from GeoServer
            projected = tempfile.NamedTemporaryFile(
                prefix="reprojected_", suffix=".json", delete=False, dir=self.workdir
            ).name
            generic_vector_reproject(
                vector_file, projected, source_crs=vec_crs, target_crs=NALCMS_PROJ4
            )

            bbox = get_bbox(projected)
            raster_url = "public:CEC_NALCMS_LandUse_2010"
            raster_bytes = geoserver.get_raster_wcs(
                bbox, geographic=False, layer=raster_url
            )
            raster_file = tempfile.NamedTemporaryFile(
                prefix="wcs_", suffix=".tiff", delete=False, dir=self.workdir
            ).name
            with open(raster_file, "wb") as f:
                f.write(raster_bytes)

        response.update_status("Accessed raster", status_percentage=20)

        categories = SIMPLE_CATEGORIES if simple_categories else TRUE_CATEGORIES
        summary_stats = SUMMARY_ZONAL_STATS

        try:
            stats = zonal_stats(
                projected,
                raster_file,
                stats=summary_stats,
                band=band,
                categorical=True,
                all_touched=touches,
                geojson_out=True,
                raster_out=False,
            )

            response.update_status("Statistic calculated", status_percentage=70)

            land_use = list()
            for stat in stats:
                lu = defaultdict(int)
                prop = stat["properties"]

                # Rename/aggregate land-use categories
                for k, v in categories.items():
                    lu[v] += prop.get(k, 0)

                prop.update(lu)
                land_use.append(lu)

            feature_collect = {"type": "FeatureCollection", "features": stats}
            response.outputs["features"].data = json.dumps(feature_collect)
            response.outputs["statistics"].data = json.dumps(land_use)

        except Exception as e:
            msg = f"Failed to perform zonal statistics using {shape_url} and {raster_url}: {e}"
            LOGGER.error(msg)
            raise Exception(msg) from e

        return response
