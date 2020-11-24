from .wps_ostrich import OstrichProcess
from raven.models import HBVEC_OST
from . import wpsio as wio
import logging
from pywps import FORMATS, LiteralInput, ComplexOutput
from pathlib import Path

LOGGER = logging.getLogger("PYWPS")

"""
Notes
-----
The configuration files for a OSTRICH calibration of the HBV-EC model and in models/ostrich-hbv-ec.
All parameters that could potentially be user-defined are tagged using {}. These tags need to be replaced by
actual values before the model is launched.
"""

params_defaults = HBVEC_OST.params(par_x01=0.05984519,
                                   par_x02=4.072232,
                                   par_x03=2.001574,
                                   par_x04=0.03473693,
                                   par_x05=0.09985144,
                                   par_x06=0.506052,
                                   par_x07=3.438486,
                                   par_x08=38.32455,
                                   par_x09=0.4606565,
                                   par_x10=0.06303738,
                                   par_x11=2.277781,
                                   par_x12=4.873686,
                                   par_x13=0.5718813,
                                   par_x14=0.04505643,
                                   par_x15=0.877607,
                                   par_x16=18.94145,
                                   par_x17=2.036937,
                                   par_x18=0.4452843,
                                   par_x19=0.6771759,
                                   par_x20=1.141608,
                                   par_x21=1.024278)

Lparams_defaults = HBVEC_OST.params(par_x01=-3.0,
                                    par_x02=0.0,
                                    par_x03=0.0,
                                    par_x04=0.0,
                                    par_x05=0.0,
                                    par_x06=0.3,
                                    par_x07=0.0,
                                    par_x08=0.0,
                                    par_x09=0.01,
                                    par_x10=0.05,
                                    par_x11=0.01,
                                    par_x12=0.0,
                                    par_x13=0.0,
                                    par_x14=0.0,
                                    par_x15=0.0,
                                    par_x16=0.0,
                                    par_x17=0.01,
                                    par_x18=0.0,
                                    par_x19=0.05,
                                    par_x20=0.8,
                                    par_x21=0.8)

Uparams_defaults = HBVEC_OST.params(par_x01=3.0,
                                    par_x02=8.0,
                                    par_x03=8.0,
                                    par_x04=0.1,
                                    par_x05=1.0,
                                    par_x06=1.0,
                                    par_x07=7.0,
                                    par_x08=100.0,
                                    par_x09=1.0,
                                    par_x10=0.1,
                                    par_x11=6.0,
                                    par_x12=5.0,
                                    par_x13=5.0,
                                    par_x14=0.2,
                                    par_x15=1.0,
                                    par_x16=30.0,
                                    par_x17=3.0,
                                    par_x18=2.0,
                                    par_x19=1.0,
                                    par_x20=1.5,
                                    par_x21=1.5)

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


class OstrichHBVECProcess(OstrichProcess):
    """
    OSTRICH emulator for the HBV-EC model.

    This process calibrates the HBV-EC model using a OSTRICH emulator. Users need to provide netCDF input
    files for
    rain, snow minimum and maximum temperature as well as potential evapotranspiration. To run diagnostics, observed
    stream flows are also required.
    """
    identifier = 'ostrich-hbv-ec'
    abstract = 'OSTRICH calibration of RAVEN HBV-EC hydrological model'
    title = ''
    version = ''
    model_cls = HBVEC_OST
    tuple_inputs = {'lowerBounds': HBVEC_OST.params,
                    'upperBounds': HBVEC_OST.params}
    inputs = [wio.ts, wio.nc_spec, wio.nc_index, lowerBounds, upperBounds, wio.algorithm, wio.max_iterations,
              wio.start_date, wio.end_date, wio.duration, wio.run_name, wio.name, wio.area,
              wio.latitude, wio.longitude, wio.elevation,
              wio.random_seed, wio.suppress_output, wio.rain_snow_fraction, wio.evaporation, wio.ow_evaporation]

    keywords = ["Ostrich", "Calibration", "DDS"]
