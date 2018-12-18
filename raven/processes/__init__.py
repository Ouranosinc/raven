from .wps_raven import RavenProcess
from .wps_gr4j_cemaneige import GR4JCemaNeigeProcess
from .wps_raven_gr4j_cemaneige import RavenGR4JCemaNeigeProcess
from .wps_raven_mohyse import RavenMOHYSEProcess
from .wps_raven_hmets import RavenHMETSProcess
from .wps_raven_hbv_ec import RavenHBVECProcess
from .wps_shape_area import ShapeAreaProcess
from .wps_shape_selection import ShapeSelectionProcess
from .wps_zonal_stats import ZonalStatisticsProcess

processes = [
    RavenProcess(),
    GR4JCemaNeigeProcess(),
    RavenGR4JCemaNeigeProcess(),
    RavenMOHYSEProcess(),
    RavenHMETSProcess(),
    RavenHBVECProcess(),
    ShapeAreaProcess(),
    ShapeSelectionProcess(),
    ZonalStatisticsProcess(),
]
