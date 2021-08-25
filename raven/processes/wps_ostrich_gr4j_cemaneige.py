import logging
from dataclasses import astuple, fields

from pywps import LiteralInput
from ravenpy.models import GR4JCN_OST

from . import wpsio as wio
from .wps_ostrich import OstrichProcess

LOGGER = logging.getLogger("PYWPS")

"""
Notes
-----
The configuration files for a OSTRICH calibration of the GR4J-Cemaneige model and in models/ostrich-gr4j-cemaneige.
All parameters that could potentially be user-defined are tagged using {}. These tags need to be replaced by
actual values before the model is launched.
"""

params_defaults = GR4JCN_OST.Params(
    GR4J_X1=0.529,
    GR4J_X2=-3.396,
    GR4J_X3=407.29,
    GR4J_X4=1.072,
    CEMANEIGE_X1=16.9,
    CEMANEIGE_X2=0.947,
)

Uparams_defaults = GR4JCN_OST.Params(
    GR4J_X1=0.9,
    GR4J_X2=0.0,
    GR4J_X3=500.0,
    GR4J_X4=1.1,
    CEMANEIGE_X1=20.0,
    CEMANEIGE_X2=1.0,
)

Lparams_defaults = GR4JCN_OST.Params(
    GR4J_X1=0.1,
    GR4J_X2=-5.0,
    GR4J_X3=100.0,
    GR4J_X4=1.0,
    CEMANEIGE_X1=10,
    CEMANEIGE_X2=0.1,
)

upperBounds = LiteralInput(
    "upperBounds",
    "Comma separated list of model parameters Upper Bounds",
    abstract="UParameters: " + ", ".join(f.name for f in fields(Uparams_defaults)),
    data_type="string",
    default=", ".join(map(str, astuple(Uparams_defaults))),
    min_occurs=0,
)

lowerBounds = LiteralInput(
    "lowerBounds",
    "Comma separated list of model parameters Lower Bounds",
    abstract="LParameters: " + ", ".join(f.name for f in fields(Lparams_defaults)),
    data_type="string",
    default=", ".join(map(str, astuple(Lparams_defaults))),
    min_occurs=0,
)


class OstrichGR4JCemaNeigeProcess(OstrichProcess):
    """
    OSTRICH emulator for the GR4J-Cemaneige model.

    This process calibrates the GR4J-Cemaneige model using a OSTRICH emulator. Users need to provide netCDF input
    files for
    rain, snow minimum and maximum temperature as well as potential evapotranspiration. To run diagnostics, observed
    stream flows are also required.
    """

    identifier = "ostrich-gr4j-cemaneige"
    abstract = "OSTRICH calibration of RAVEN GR4J + CEMANEIGE hydrological model"
    title = ""
    version = ""
    model_cls = GR4JCN_OST
    tuple_inputs = {
        "lowerBounds": GR4JCN_OST.Params,
        "upperBounds": GR4JCN_OST.Params,
    }
    inputs = [
        wio.ts,
        wio.nc_spec,
        wio.nc_index,
        lowerBounds,
        upperBounds,
        wio.algorithm,
        wio.max_iterations,
        wio.start_date,
        wio.end_date,
        wio.duration,
        wio.run_name,
        wio.area,
        wio.latitude,
        wio.longitude,
        wio.elevation,
        wio.random_seed,
        wio.random_numbers,
        wio.suppress_output,
        wio.evaporation,
        wio.rain_snow_fraction,
    ]

    keywords = ["Ostrich", "Calibration", "DDS"]
