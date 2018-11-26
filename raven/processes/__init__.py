from .wps_wordcounter import WordCounter
from .wps_inout import InOut
from .wps_sleep import Sleep
from .wps_gr4j_cemaneige import GR4JCemaNeigeProcess
# from .wps_raven_gr4j_cemaneige import RavenGR4JCemaNeigeProcess
# from .wps_raven_hmets import RavenHMETSProcess

processes = [
    GR4JCemaNeigeProcess(),
    # RavenGR4JCemaNeigeProcess(),
    WordCounter(),
    InOut(),
    Sleep(),
]
