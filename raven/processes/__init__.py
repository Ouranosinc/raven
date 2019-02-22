from .wps_raven import RavenProcess
from .wps_gr4j_cemaneige import GR4JCemaNeigeProcess
from .wps_raven_gr4j_cemaneige import RavenGR4JCemaNeigeProcess
from .wps_raven_mohyse import RavenMOHYSEProcess
from .wps_raven_hmets import RavenHMETSProcess
from .wps_raven_hbv_ec import RavenHBVECProcess
from .wps_generic_shape_properties import ShapeAreaProcess
from .wps_hydrobasins_shape_selection import ShapeSelectionProcess
from .wps_generic_zonal_stats import ZonalStatisticsProcess
from .wps_generic_raster_subset import RasterSubsetProcess
from .wps_ostrich_gr4j_cemaneige import OstrichGR4JCemaNeigeProcess
from .wps_ostrich_mohyse import OstrichMOHYSEProcess
from .wps_ostrich_hmets import OstrichHMETSProcess
from .wps_ostrich_hbv_ec import OstrichHBVECProcess
from .wps_objective_functions import ObjectiveFunctionProcess
from .wps_regionalisation import RegionalisationProcess


processes = [
    RavenProcess(),
    GR4JCemaNeigeProcess(),
    RasterSubsetProcess(),
    RavenGR4JCemaNeigeProcess(),
    RavenMOHYSEProcess(),
    RavenHMETSProcess(),
    RavenHBVECProcess(),
    ShapeAreaProcess(),
    ShapeSelectionProcess(),
    ZonalStatisticsProcess(),
    OstrichGR4JCemaNeigeProcess(),
    OstrichMOHYSEProcess(),
    OstrichHMETSProcess(),
    OstrichHBVECProcess(),
    ObjectiveFunctionProcess(),
    # RegionalisationProcess(),
]
