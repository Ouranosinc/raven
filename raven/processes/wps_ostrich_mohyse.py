from .wps_ostrich import OstrichProcess
from raven.models import MOHYSE_OST
from . import wpsio as wio
import logging
from pywps import FORMATS, LiteralInput, ComplexOutput
from pathlib import Path

LOGGER = logging.getLogger("PYWPS")

"""
Notes
-----
The configuration files for a OSTRICH calibration of the MOHYSE model and in models/ostrich-mohyse.
All parameters that could potentially be user-defined are tagged using {}. These tags need to be replaced by
actual values before the model is launched.
"""

params_defaults = MOHYSE_OST.params(par_x01=1.0000,
                                    par_x02=0.0468,
                                    par_x03=4.2952,
                                    par_x04=2.6580,
                                    par_x05=0.4038,
                                    par_x06=0.0621,
                                    par_x07=0.0273,
                                    par_x08=0.0453,
                                    par_x09=0.9039,
                                    par_x10=5.6167)

Lparams_defaults = MOHYSE_OST.params(par_x01=0.01,
                                     par_x02=0.01,
                                     par_x03=0.01,
                                     par_x04=-5.00,
                                     par_x05=0.01,
                                     par_x06=0.01,
                                     par_x07=0.01,
                                     par_x08=0.01,
                                     par_x09=0.01,
                                     par_x10=0.01
                                     )
Uparams_defaults = MOHYSE_OST.params(par_x01=20.0,
                                     par_x02=1.0,
                                     par_x03=20.0,
                                     par_x04=5.0,
                                     par_x05=0.5,
                                     par_x06=1.0,
                                     par_x07=1.0,
                                     par_x08=1.0,
                                     par_x09=15.0,
                                     par_x10=15.0
                                     )

upperBounds = LiteralInput('upperBounds', 'Comma separated list of model parameters Upper Bounds',
                           abstract='UParameters: ' + ', '.join(Uparams_defaults._fields),
                           data_type='string',
                           default=', '.join(str(p) for p in list(Uparams_defaults)),
                           min_occurs=0)

lowerBounds = LiteralInput('lowerBounds', 'Comma separated list of model parameters Lower Bounds',
                           abstract='LParameters: ' + ', '.join(Lparams_defaults._fields),
                           data_type='string',
                           default=', '.join(str(p) for p in list(Lparams_defaults)),
                           min_occurs=0)


class OstrichMOHYSEProcess(OstrichProcess):
    """
    OSTRICH emulator for the MOHYSE model.

    This process calibrates the MOHYSE model using a OSTRICH emulator. Users need to provide netCDF input
    files for
    rain, snow minimum and maximum temperature as well as potential evapotranspiration. To run diagnostics, observed
    stream flows are also required.
    """
    identifier = 'ostrich-mohyse'
    abstract = 'OSTRICH calibration of RAVEN MOHYSE hydrological model'
    title = ''
    version = ''
    model_cls = MOHYSE_OST
    tuple_inputs = {'lowerBounds': MOHYSE_OST.params,
                    'upperBounds': MOHYSE_OST.params}
    inputs = [wio.ts, wio.nc_spec, lowerBounds, upperBounds, wio.algorithm,
              wio.max_iterations, wio.start_date, wio.end_date,
              wio.duration, wio.run_name, wio.name, wio.area, wio.latitude, wio.longitude, wio.elevation,
              wio.random_seed, wio.suppress_output, wio.rain_snow_fraction, wio.evaporation]

    keywords = ["Ostrich", "Calibration", "DDS"]
