"""Processes for the Raven WPS server."""

from .wps_generic_raster_subset import RasterSubsetProcess
from .wps_generic_shape_properties import ShapePropertiesProcess
from .wps_generic_terrain_analysis import TerrainAnalysisProcess
from .wps_generic_zonal_stats import ZonalStatisticsProcess
from .wps_hydrobasins_shape_selection import HydroBasinsSelectionProcess
from .wps_nalcms_zonal_stats import NALCMSZonalStatisticsProcess
from .wps_nalcms_zonal_stats_raster import NALCMSZonalStatisticsRasterProcess

geo = [
    RasterSubsetProcess(),
    ShapePropertiesProcess(),
    HydroBasinsSelectionProcess(),
    TerrainAnalysisProcess(),
    ZonalStatisticsProcess(),
    NALCMSZonalStatisticsProcess(),
    NALCMSZonalStatisticsRasterProcess(),
]


processes = geo
