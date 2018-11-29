from .wps_wordcounter import WordCounter
from .wps_inout import InOut
from .wps_sleep import Sleep
from .wps_raven import RavenProcess
from .wps_gr4j_cemaneige import GR4JCemaNeigeProcess
from .wps_raven_gr4j_cemaneige import RavenGR4JCemaNeigeProcess
from .wps_raven_hmets import RavenHMETSProcess

processes = [
    RavenProcess(),
    GR4JCemaNeigeProcess(),
    RavenGR4JCemaNeigeProcess(),
    RavenHMETSProcess(),
    WordCounter(),
    InOut(),
    Sleep(),
]
