import datetime as dt
import os
from .wps_raven import RavenProcess
from pywps import LiteralInput
from raven.models import GR4JCemaneige
from . import wpsio as wio

import logging
from collections import OrderedDict as Odict

LOGGER = logging.getLogger("PYWPS")


"""
Notes
-----

The configuration files for RAVEN's GR4J-Cemaneige model and in models/raven-gr4j-cemaneige.
All parameters that could potentially be user-defined are tagged using {}. These tags need to be replaced by
actual values before the model is launched. 
"""

param_defaults = Odict([('SOIL_PROD', 0.696),
                   ('GR4J_X2', 0.7),
                   ('GR4J_X3', 19.7),
                   ('GR4J_X4', 2.09),
                   ('AvgAnnualSnow', 123.3),
                   ('AirSnowCoeff', 0.75)])

params = LiteralInput('params', 'Comma separated list of model parameters',
                      abstract='Parameters: ' + ', '.join(param_defaults.keys()),
                      data_type='string',
                      default=', '.join(str(p) for p in param_defaults.values()),
                      min_occurs=0)

init = LiteralInput('init', 'Initial soil conditions',
                    abstract='Underground reservoir levels: SOIL_0, SOIL_1',
                    data_type='string',
                    default='0, 0',
                    min_occurs=0)

class RavenGR4JCemaNeigeProcess(RavenProcess):
    """
    RAVEN emulator for the GR4J-Cemaneige model.

    This process runs the GR4J-Cemaneige model using a RAVEN emulator. Users need to provide netCDF input files for
    rain, snow minimum and maximum temperature as well as potential evapotranspiration. To run diagnostics, observed
    stream flows are also required.


    """
    identifier = 'raven-gr4j-cemaneige'
    abstract = 'Raven GR4J + CEMANEIGE hydrological model'
    title = ''
    version = ''
    model_cls = GR4JCemaneige
    param_arrays = ['params', 'init']

    inputs = [wio.ts, params, wio.start_date, wio.end_date, wio.duration, init, wio.run_name,
              wio.name, wio.area, wio.latitude, wio.longitude, wio.elevation]

