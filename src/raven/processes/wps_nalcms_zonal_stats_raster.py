import json
import logging
import tempfile
from collections import defaultdict
from pathlib import Path

import geopandas as gpd
from pywps import FORMATS, ComplexOutput, Process
from pywps.inout.outputs import MetaFile, MetaLink4
from rasterstats import zonal_stats

from raven.utilities.checks import single_file_check
from raven.utilities.geo import generic_vector_reproject
from raven.utilities.io import archive_sniffer, crs_sniffer, raster_datatype_sniffer
from raven.utils import (
    NALCMS_PROJ4,
    SIMPLE_CATEGORIES,
    SUMMARY_ZONAL_STATS,
    TRUE_CATEGORIES,
    gather_dem_tile,
    zonalstats_raster_file,
)

from . import wpsio as wio

LOGGER = logging.getLogger("PYWPS")


class NALCMSZonalStatisticsRasterProcess(Process):
    """Given files containing vector data and raster data, perform zonal statistics of the overlapping regions"""

    def __init__(self):
        inputs = [
            wio.shape,
            wio.land_use_raster,
            wio.simple_categories,
            wio.raster_band,
            wio.select_all_touching,
        ]

        outputs = [
            wio.features,
            wio.statistics,
            ComplexOutput(
                "raster",
                "DEM grid subset by the requested shape.",
                abstract="Zipped raster grid(s) of land-use using either standard or simplified UNFAO categories.",
                as_reference=True,
                supported_formats=[FORMATS.META4],
            ),
        ]

        super().__init__(
            self._handler,
            identifier="nalcms-zonal-stats-raster",
            title="NALCMS Land Use Zonal Statistics with raster output",
            version="1.0",
            abstract="Return zonal statistics, land-use cover, and raster grid for the CEC NALCMS based on the boundaries of a vector file.",
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

        # For raster files using the UNFAO Land Cover Classification System (19 types)
        if "raster" in request.inputs:
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

        else:
            raster_url = None
            # using the NALCMS data from GeoServer
            projected = tempfile.NamedTemporaryFile(
                prefix="reprojected_", suffix=".json", delete=False, dir=self.workdir
            ).name
            generic_vector_reproject(
                vector_file, projected, source_crs=vec_crs, target_crs=NALCMS_PROJ4
            )
            gdf = gpd.read_file(projected)
            if sum(gdf.area) / 1e6 > 1e5:
                LOGGER.warning(f"Vector shape has area of {sum(gdf.area) / 1e6} km2.")
                raise Exception(
                    "NALCMS zonal statistics only supported for areas smaller than 100,000 km2."
                )

            raster_file = gather_dem_tile(
                projected,
                self.workdir,
                geographic=False,
                raster="public:CEC_NALCMS_LandUse_2010",
            )

        data_type = raster_datatype_sniffer(raster_file)
        response.update_status("Accessed raster", status_percentage=10)

        categories = SIMPLE_CATEGORIES if simple_categories else TRUE_CATEGORIES
        summary_stats = SUMMARY_ZONAL_STATS

        try:
            # Use zonalstats to produce a GeoJSON
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

            land_use = []
            for stat in stats:
                lu = defaultdict(int)
                prop = stat["properties"]

                # Rename/aggregate land-use categories
                for k, v in categories.items():
                    # Fiona v1.9 API changes; Access to a protected method of class instance - Needs rewrite
                    lu[v] += prop.get(k, 0)

                prop.update(lu)
                land_use.append(lu)
                # prop['mini_raster_array'] = pickle.dumps(prop['mini_raster_array'], protocol=0).decode()

            # Use zonalstats to produce sets of raster grids
            raster_subset = zonal_stats(
                projected,
                raster_file,
                stats=summary_stats,
                band=band,
                categorical=True,
                all_touched=touches,
                geojson_out=False,
                raster_out=True,
            )

            raster_out = zonalstats_raster_file(
                raster_subset,
                working_dir=self.workdir,
                data_type=data_type,
                crs=NALCMS_PROJ4,
                zip_archive=False,
            )

            ml = MetaLink4(
                "rasters_out",
                "Metalink to series of GeoTIFF raster files",
                workdir=self.workdir,
            )
            for raster in raster_out:
                mf = MetaFile(Path(raster).name, "Raster subset", fmt=FORMATS.GEOTIFF)
                mf.file = raster
                ml.append(mf)

            feature_collect = {"type": "FeatureCollection", "features": stats}
            response.outputs["features"].data = json.dumps(feature_collect)
            response.outputs["statistics"].data = json.dumps(land_use)
            response.outputs["raster"].data = ml.xml

        except Exception as e:
            msg = f"Failed to perform raster subset using {shape_url}{f' and {raster_url} ' if raster_url else ''}: {e}"
            LOGGER.error(msg)
            raise Exception(msg) from e

        return response
