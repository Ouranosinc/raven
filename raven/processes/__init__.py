from .wps_raven import RavenProcess
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
from .wps_nalcms_zonal_stats_raster import NALCMSZonalStatisticsRasterProcess
from .wps_ostrich_gr4j_cemaneige import OstrichGR4JCemaNeigeProcess
from .wps_ostrich_mohyse import OstrichMOHYSEProcess
from .wps_ostrich_hmets import OstrichHMETSProcess
from .wps_ostrich_hbv_ec import OstrichHBVECProcess
from .wps_objective_functions import ObjectiveFunctionProcess
from .wps_regionalisation import RegionalisationProcess
from .wps_raven_multi_model import RavenMultiModelProcess
from .wps_graph_ensemble_uncertainty import GraphEnsUncertaintyProcess
from .wps_graph_single_hydrograph import GraphSingleHydrographProcess
from .wps_indicator_analysis import GraphIndicatorAnalysis
from .wps_graph_objective_function_fit import GraphObjectiveFunctionFitProcess
from .wps_graph_fit import GraphFitProcess
from .wps_forecast_flood_risk import ForecastFloodRiskProcess
from .wps_climatology_esp import ClimatologyEspProcess
from .wps_forecast_evaluation import HindcastEvaluationProcess
from .wps_hindcast import HindcastingProcess
from .wps_realtime_forecast import RealtimeForecastProcess
from .wps_graph_forecast_uncertainty import GraphFcstUncertaintyProcess


modeling = [
    RavenProcess(),
    RavenGR4JCemaNeigeProcess(),
    RavenMOHYSEProcess(),
    RavenHMETSProcess(),
    RavenHBVECProcess(),
    RavenMultiModelProcess(),
    ClimatologyEspProcess(),
    HindcastEvaluationProcess(),
    HindcastingProcess(),
    RealtimeForecastProcess(),
]

calibration = [
    OstrichGR4JCemaNeigeProcess(),
    OstrichMOHYSEProcess(),
    OstrichHMETSProcess(),
    OstrichHBVECProcess(),
]

geo = [
    RasterSubsetProcess(),
    ShapePropertiesProcess(),
    HydroBasinsSelectionProcess(),
    TerrainAnalysisProcess(),
    ZonalStatisticsProcess(),
    NALCMSZonalStatisticsProcess(),
    NALCMSZonalStatisticsRasterProcess(),
]

analysis = [
    ObjectiveFunctionProcess(),
    RegionalisationProcess(),
    ForecastFloodRiskProcess(),
]

graphs = [
    GraphEnsUncertaintyProcess(),
    GraphSingleHydrographProcess(),
    GraphIndicatorAnalysis(),
    GraphObjectiveFunctionFitProcess(),
    GraphFitProcess(),
    GraphFcstUncertaintyProcess(),
]

processes = modeling + calibration + geo + analysis + graphs
