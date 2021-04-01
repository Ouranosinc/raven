import logging

from ravenpy.models import Ostrich

from raven.processes import RavenProcess

from . import wpsio as wio

LOGGER = logging.getLogger("PYWPS")


class OstrichProcess(RavenProcess):
    identifier = "ostrich"
    abstract = "Ostrich calibration framework"
    title = (
        "Run an Ostrich calibration of a hydrologic model setup using the Raven hydrological framework with model"
        "configuration files and forcing time series. In "
        "the `rvt` file, only provide the name of the forcing file, not an absolute or relative path."
    )
    version = "0.1"
    model_cls = Ostrich
    inputs = [wio.ts, wio.conf]
    outputs = [
        wio.calibration,
        wio.hydrograph,
        wio.storage,
        wio.solution,
        wio.diagnostics,
        wio.calibparams,
        wio.rv_config,
    ]
