import json
import logging
import tempfile
from collections import defaultdict

from pywps import ComplexOutput
from pywps import LiteralInput, ComplexInput
from pywps import Process, FORMATS
from pywps.app.Common import Metadata
from rasterstats import zonal_stats

from raven.utilities import gis
from raven.utils import (
    archive_sniffer,
    crs_sniffer,
    generic_vector_reproject,
    raster_datatype_sniffer,
    single_file_check,
    zonalstats_raster_file,
)

LOGGER = logging.getLogger("PYWPS")

TRUE_CATEGORIES = {
    0: "Ocean",
    1: "Temperate or sub-polar needleleaf forest",
    2: "Sub-polar taiga needleleaf forest",
    3: "Tropical or sub-tropical broadleaf evergreen forest",
    4: "Tropical or sub-tropical broadleaf deciduous forest",
    5: "Temperate or sub-polar broadleaf deciduous forest",
    6: "Mixed forest",
    7: "Tropical or sub-tropical shrubland",
    8: "Temperate or sub-polar shrubland",
    9: "Tropical or sub-tropical grassland",
    10: "Temperate or sub-polar grassland",
    11: "Sub-polar or polar shrubland-lichen-moss",
    12: "Sub-polar or polar grassland-lichen-moss",
    13: "Sub-polar or polar barren-lichen-moss",
    14: "Wetland",
    15: "Cropland",
    16: "Barren lands",
    17: "Urban",
    18: "Water",
    19: "Snow and Ice",
}

simplified = {
    "Ocean": [0],
    "Forest": [1, 2, 3, 4, 5, 6],
    "Shrubs": [7, 8, 11],
    "Grass": [9, 10, 12, 13, 16],
    "Wetland": [14],
    "Crops": [15],
    "Urban": [17],
    "Water": [18],
    "SnowIce": [19],
}
SIMPLE_CATEGORIES = {i: cat for (cat, ids) in simplified.items() for i in ids}

SUMMARY_ZONAL_STATS = ["count", "nodata", "nan"]
NALCMS_PROJ4 = (
    "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs=True"
)


class NALCMSZonalStatisticsRasterProcess(Process):
    """Given files containing vector data and raster data, perform zonal statistics of the overlapping regions"""

    def __init__(self):
        inputs = [
            ComplexInput(
                "shape",
                "Vector Shape",
                abstract="An ESRI Shapefile, GML, JSON, GeoJSON, or single layer GeoPackage."
                " The ESRI Shapefile must be zipped and contain the .shp, .shx, and .dbf.",
                min_occurs=1,
                max_occurs=1,
                supported_formats=[
                    FORMATS.GEOJSON,
                    FORMATS.GML,
                    FORMATS.JSON,
                    FORMATS.SHP,
                ],
            ),
            ComplexInput(
                "raster",
                "Gridded Land Use raster data set",
                abstract="The Land Use raster to be queried. Default is the CEC NALCMS 2010. Provided raster "
                "must use the UN FAO Land Cover Classification System (19 types).",
                metadata=[
                    Metadata(
                        "Commission for Environmental Cooperation North American Land Change Monitoring System",
                        "http://www.cec.org/tools-and-resources/map-files/land-cover-2010-landsat-30m",
                    ),
                    Metadata(
                        "Latifovic, R., Homer, C., Ressl, R., Pouliot, D., Hossain, S.N., Colditz, R.R.,"
                        "Olthof, I., Giri, C., Victoria, A., (2012). North American land change "
                        "monitoring system. In: Giri, C., (Ed), Remote Sensing of Land Use and Land "
                        "Cover: Principles and Applications, CRC-Press, pp. 303-324"
                    ),
                ],
                min_occurs=0,
                max_occurs=1,
                supported_formats=[FORMATS.GEOTIFF],
            ),
            LiteralInput(
                "simple_categories",
                "Use simplified land classification categories for hydrological "
                "modeling purposes.",
                data_type="boolean",
                default="false",
                min_occurs=0,
                max_occurs=1,
            ),
            LiteralInput(
                "band",
                "Raster band",
                data_type="integer",
                default=1,
                abstract="Band of raster examined to perform zonal statistics.",
                min_occurs=0,
                max_occurs=1,
            ),
            LiteralInput(
                "select_all_touching",
                "Additionally select boundary pixels that are touched by shape.",
                data_type="boolean",
                default="false",
                min_occurs=0,
                max_occurs=1,
            ),
        ]

        outputs = [
            ComplexOutput(
                "features",
                "DEM properties within the region defined by the vector provided.",
                abstract="Category pixel counts using either standard or simplified UNFAO categories",
                supported_formats=[FORMATS.GEOJSON],
            ),
            ComplexOutput(
                "statistics",
                "DEM properties by feature",
                abstract="Land-use type pixel counts using either standard or simplified UNFAO categories.",
                supported_formats=[FORMATS.JSON],
            ),
            ComplexOutput(
                "raster",
                "DEM grid subset by the requested shape.",
                abstract="Zipped raster grid(s) of land-use using either standard or simplified UNFAO categories.",
                supported_formats=[FORMATS.ZIP],
            ),
        ]

        super(NALCMSZonalStatisticsRasterProcess, self).__init__(
            self._handler,
            identifier="nalcms-zonal-stats-raster",
            title="NALCMS Land Use Zonal Statistics with raster output",
            version="1.0",
            abstract="Return zonal statistics, land-use cover, and raster grid for the CEC NALCMS based "
            "on the boundaries of a vector file.",
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
                msg = (
                    f"CRS for files {vector_file} and {raster_file} are not the same. Reprojecting..."
                )
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

            bbox = gis.get_bbox(projected)
            raster_url = "public:CEC_NALCMS_LandUse_2010"
            raster_bytes = gis.get_raster_wcs(bbox, geographic=False, layer=raster_url)
            raster_file = tempfile.NamedTemporaryFile(
                prefix="wcs_", suffix=".tiff", delete=False, dir=self.workdir
            ).name
            with open(raster_file, "wb") as f:
                f.write(raster_bytes)

        data_type = raster_datatype_sniffer(raster_file)
        response.update_status("Accessed raster", status_percentage=10)

        if simple_categories:
            categories = SIMPLE_CATEGORIES
        else:
            categories = TRUE_CATEGORIES
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

            land_use = list()
            for stat in stats:
                lu = defaultdict(lambda: 0)
                prop = stat["properties"]

                # Rename/aggregate land-use categories
                for k, v in categories.items():
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
                crs=vec_crs,
            )

            feature_collect = {"type": "FeatureCollection", "features": stats}
            response.outputs["features"].data = json.dumps(feature_collect)
            response.outputs["statistics"].data = json.dumps(land_use)
            response.outputs["raster"].data = raster_out

        except Exception as e:
            msg = f"Failed to perform zonal statistics using {shape_url} and {raster_url}: {e}"
            LOGGER.error(msg)
            raise Exception(msg) from e

        return response
