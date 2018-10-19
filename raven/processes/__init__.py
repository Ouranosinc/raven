from .wps_wordcounter import WordCounter
from .wps_inout import InOut
from .wps_sleep import Sleep
from .wps_gr4j_cemaneige import GR4JCemaNeigeProcess

processes = [
    GR4JCemaNeigeProcess(),
    WordCounter(),
    InOut(),
    Sleep(),
]
