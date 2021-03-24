from raven import config
from .wps_raven import RavenProcess
from ravenpy.models import GR4JCN
from . import wpsio as wio
import logging
from pywps import LiteralInput

LOGGER = logging.getLogger("PYWPS")

"""
Notes
-----

The configuration files for RAVEN's GR4J-Cemaneige model and in models/raven-gr4j-cemaneige.
All parameters that could potentially be user-defined are tagged using {}. These tags need to be replaced by
actual values before the model is launched.
"""

params_defaults = GR4JCN.params(
    GR4J_X1=0.529,
    GR4J_X2=-3.396,
    GR4J_X3=407.29,
    GR4J_X4=1.072,
    CEMANEIGE_X1=16.9,
    CEMANEIGE_X2=0.947,
)

params = LiteralInput(
    "params",
    "Comma separated list of model parameters",
    abstract="Parameters: " + ", ".join(params_defaults._fields),
    data_type="string",
    default=", ".join(str(p) for p in list(params_defaults)),
    min_occurs=0,
    max_occurs=config.max_parallel_processes,
)


class RavenGR4JCemaNeigeProcess(RavenProcess):
    """
    RAVEN emulator for the GR4J-Cemaneige model.

    This process runs the GR4J-Cemaneige model using a RAVEN emulator. Users need to provide netCDF input files for
    rain, snow minimum and maximum temperature as well as potential evapotranspiration. To run diagnostics, observed
    stream flows are also required.


    """

    identifier = "raven-gr4j-cemaneige"
    abstract = "Raven GR4J + CEMANEIGE hydrological model"
    title = ""
    version = ""
    model_cls = GR4JCN
    tuple_inputs = {"params": GR4JCN.params}

    inputs = [
        wio.ts,
        wio.nc_spec,
        params,
        wio.start_date,
        wio.end_date,
        wio.nc_index,
        wio.duration,
        wio.run_name,
        wio.name,
        wio.area,
        wio.latitude,
        wio.longitude,
        wio.elevation,
        wio.evaporation,
        wio.rain_snow_fraction,
        wio.rvc,
    ]
