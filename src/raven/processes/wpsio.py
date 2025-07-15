"""This module contains the WPS inputs and outputs that are reused across multiple WPS processes."""

from pywps import FORMATS, ComplexInput, ComplexOutput, LiteralInput
from pywps.app.Common import Metadata

# ---------------------------------------- #
# ---------------- Inputs ---------------- #
# ---------------------------------------- #

# --- GIS Inputs --- #

region_vector = ComplexInput(
    "region_vector",
    "Vector shape file of a region",
    abstract="An ESRI Shapefile, GML, JSON, GeoJSON, or single layer GeoPackage."
    " The ESRI Shapefile must be zipped and contain the .shp, .shx, and .dbf.",
    min_occurs=1,
    max_occurs=1,
    supported_formats=[
        FORMATS.GEOJSON,
        FORMATS.GML,
        FORMATS.JSON,
        FORMATS.SHP,
        FORMATS.ZIP,
    ],
)

shape = ComplexInput(
    "shape",
    "Vector shape of a region",
    abstract="An ESRI Shapefile, GML, JSON, GeoJSON, or single layer GeoPackage."
    " The ESRI Shapefile must be zipped and contain the .shp, .shx, and .dbf.",
    min_occurs=1,
    max_occurs=1,
    supported_formats=[
        FORMATS.GEOJSON,
        FORMATS.GML,
        FORMATS.JSON,
        FORMATS.SHP,
        FORMATS.ZIP,
    ],
)

land_use_raster = ComplexInput(
    "raster",
    "Gridded Land Use raster data set",
    abstract="The Land Use raster to be queried. Default is the CEC NALCMS 2010. Provided "
    "raster "
    "must use the UN FAO Land Cover Classification System (19 types).",
    metadata=[
        Metadata(
            "Commission for Environmental Cooperation North American Land Change Monitoring System",
            "https://www.cec.org/tools-and-resources/map-files/land-cover-2010-landsat-30m",
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
)

dem_raster = ComplexInput(
    "raster",
    "Gridded raster data set",
    abstract="The DEM to be queried. Defaults to the EarthEnv-DEM90 product.",
    metadata=[
        Metadata("EarthEnv-DEM90", "https://www.earthenv.org/DEM"),
        Metadata(
            "Robinson, Natalie, James Regetz, and Robert P. Guralnick (2014). "
            "EarthEnv-DEM90: A Nearly-Global, Void-Free, Multi-Scale Smoothed, 90m Digital "
            "Elevation Model from Fused ASTER and SRTM Data. ISPRS Journal of "
            "Photogrammetry and Remote Sensing 87: 57–67.",
            "https://doi.org/10.1016/j.isprsjprs.2013.11.002",
        ),
    ],
    min_occurs=0,
    max_occurs=1,
    supported_formats=[FORMATS.GEOTIFF],
)

simple_categories = LiteralInput(
    "simple_categories",
    "Use simplified land classification categories for hydrological modeling purposes.",
    data_type="boolean",
    default="false",
    min_occurs=0,
    max_occurs=1,
)

raster_band = LiteralInput(
    "band",
    "Raster band",
    data_type="integer",
    default=1,
    abstract="Band of raster examined to perform zonal statistics.",
    min_occurs=0,
    max_occurs=1,
)

select_all_touching = LiteralInput(
    "select_all_touching",
    "Additionally select boundary pixels that are touched by shape.",
    data_type="boolean",
    default="false",
    min_occurs=0,
    max_occurs=1,
)

features = ComplexOutput(
    "features",
    "DEM properties within the region defined by the vector provided.",
    abstract="Category pixel counts using either standard or simplified UNFAO categories",
    supported_formats=[FORMATS.GEOJSON],
)
statistics = ComplexOutput(
    "statistics",
    "DEM properties by feature",
    abstract="Land-use type pixel counts using either standard or simplified UNFAO categories.",
    supported_formats=[FORMATS.JSON],
)


# TODO: Add configuration files to output
# config = ComplexOutput('config', 'Configuration files',
#                        abstract="Link to configuration files.",
#                        supported_formats=)
