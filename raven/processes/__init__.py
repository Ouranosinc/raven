from .wps_raven import RavenProcess
from .wps_gr4j_cemaneige import GR4JCemaNeigeProcess
from .wps_raven_gr4j_cemaneige import RavenGR4JCemaNeigeProcess
from .wps_raven_mohyse import RavenMOHYSEProcess
from .wps_raven_hmets import RavenHMETSProcess
from .wps_raven_hbv_ec import RavenHBVECProcess
from .wps_generic_shape_properties import ShapePropertiesProcess
from .wps_hydrobasins_shape_selection import HydroBasinsSelectionProcess
from .wps_generic_raster_subset import RasterSubsetProcess
from .wps_generic_terrain_analysis import TerrainAnalysisProcess
from .wps_generic_zonal_stats import ZonalStatisticsProcess
from .wps_nalcms_zonal_stats import NALCMSZonalStatisticsProcess
from .wps_ostrich_gr4j_cemaneige import OstrichGR4JCemaNeigeProcess
from .wps_ostrich_mohyse import OstrichMOHYSEProcess
from .wps_ostrich_hmets import OstrichHMETSProcess
from .wps_ostrich_hbv_ec import OstrichHBVECProcess
from .wps_objective_functions import ObjectiveFunctionProcess
from .wps_regionalisation import RegionalisationProcess
from .wps_raven_multi_model import RavenMultiModelProcess
from .wps_graph_ensemble_uncertainty import GraphEnsUncertaintyProcess
from .wps_graph_single_hydrograph import GraphSingleHydrographProcess
from .wps_q_stats import TSStatsProcess, FreqAnalysisProcess, FitProcess, BaseFlowIndexProcess
from .wps_indicator_analysis import GraphIndicatorAnalysis
from .wps_graph_objective_function_fit import GraphObjectiveFunctionFitProcess
from .wps_graph_fit import GraphFitProcess
from .wps_forecast_flood_risk import ForecastFloodRiskProcess

processes = [
    RavenProcess(),
    GR4JCemaNeigeProcess(),
    RasterSubsetProcess(),
    RavenGR4JCemaNeigeProcess(),
    RavenMOHYSEProcess(),
    RavenHMETSProcess(),
    RavenHBVECProcess(),
    ShapePropertiesProcess(),
    HydroBasinsSelectionProcess(),
    TerrainAnalysisProcess(),
    ZonalStatisticsProcess(),
    NALCMSZonalStatisticsProcess(),
    OstrichGR4JCemaNeigeProcess(),
    OstrichMOHYSEProcess(),
    OstrichHMETSProcess(),
    OstrichHBVECProcess(),
    ObjectiveFunctionProcess(),
    RavenMultiModelProcess(),
    GraphEnsUncertaintyProcess(),
    GraphSingleHydrographProcess(),
    TSStatsProcess(),
    FitProcess(),
    FreqAnalysisProcess(),
    BaseFlowIndexProcess(),
    GraphIndicatorAnalysis(),
    RegionalisationProcess(),
    GraphObjectiveFunctionFitProcess(),
    GraphFitProcess(),
    ForecastFloodRiskProcess(),
]
