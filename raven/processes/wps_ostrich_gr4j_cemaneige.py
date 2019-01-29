from .wps_ostrich import OstrichProcess
from raven.models import GR4JCN_OST
from . import wpsio as wio
import logging
from pywps import LiteralInput

LOGGER = logging.getLogger("PYWPS")

"""
Notes
-----

The configuration files for a OSTRICH calibration of the GR4J-Cemaneige model and in models/ostrich-gr4j-cemaneige.
All parameters that could potentially be user-defined are tagged using {}. These tags need to be replaced by
actual values before the model is launched.
"""

params_defaults = GR4JCN_OST.RVP.params(GR4J_X1=0.529,
                                    GR4J_X2=-3.396,
                                    GR4J_X3=407.29,
                                    GR4J_X4=1.072,
                                    CEMANEIGE_X1=16.9,
                                    CEMANEIGE_X2=0.947)

params = LiteralInput('params', 'Comma separated list of model parameters',
                      abstract='Parameters: ' + ', '.join(params_defaults._fields),
                      data_type='string',
                      default=', '.join(str(p) for p in list(params_defaults)),
                      min_occurs=0)


class OstrichGR4JCemaNeigeProcess(OstrichProcess):
    """
    OSTRICH emulator for the GR4J-Cemaneige model.

    This process calibrates the GR4J-Cemaneige model using a OSTRICH emulator. Users need to provide netCDF input files for
    rain, snow minimum and maximum temperature as well as potential evapotranspiration. To run diagnostics, observed
    stream flows are also required.
    """

    identifier = 'ostrich-gr4j-cemaneige'
    abstract = 'OSTRICH calibration of RAVEN GR4J + CEMANEIGE hydrological model'
    title = ''
    version = ''
    model_cls = GR4JCN_OST
    tuple_inputs = {'params': GR4JCN_OST.RVP.params}

    inputs = [wio.ts, params, wio.start_date, wio.end_date, wio.duration, wio.run_name,
              wio.name, wio.area, wio.latitude, wio.longitude, wio.elevation]
