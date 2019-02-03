from .wps_ostrich import OstrichProcess
from raven.models import MOHYSE_OST
from . import wpsio as wio
import logging
from pywps import FORMATS, LiteralInput, ComplexOutput
import pdb
from pathlib import Path

LOGGER = logging.getLogger("PYWPS")

"""
Notes
-----
The configuration files for a OSTRICH calibration of the MOHYSE model and in models/ostrich-mohyse.
All parameters that could potentially be user-defined are tagged using {}. These tags need to be replaced by
actual values before the model is launched.
"""

# params_defaults = MOHYSE_OST.params(par_x01=1.0000,
#                                     par_x02=0.0468,
#                                     par_x03=4.2952,
#                                     par_x04=2.6580,
#                                     par_x05=0.4038,
#                                     par_x06=0.0621,
#                                     par_x07=0.0273,
#                                     par_x08=0.0453,
#                                     par_x09=0.9039,
#                                     par_x10=5.6167)

# Lparams_defaults = MOHYSE_OST.params(par_x01=0.01,
#                                      par_x02=0.01,
#                                      par_x03=0.01,
#                                      par_x04=-5.00,
#                                      par_x05=0.01,
#                                      par_x06=0.01,
#                                      par_x07=0.01,
#                                      par_x08=0.01,
#                                      par_x09=0.01,
#                                      par_x10=0.01)

# Uparams_defaults = MOHYSE_OST.params(par_x01=20.0,
#                                      par_x02=1.0,
#                                      par_x03=20.0,
#                                      par_x04=5.0,
#                                      par_x05=0.5,
#                                      par_x06=1.0,
#                                      par_x07=1.0,
#                                      par_x08=1.0,
#                                      par_x09=15.0,
#                                      par_x10=15.0)

params_defaults = MOHYSE_OST.params(par_x01=1.0000,
                                    par_x02=0.0468,
                                    par_x03=4.2952,
                                    par_x04=2.6580,
                                    par_x05=0.4038,
                                    par_x06=0.0621,
                                    par_x07=0.0273,
                                    par_x08=0.0453)
Lparams_defaults = MOHYSE_OST.params(par_x01=0.01,
                                     par_x02=0.01,
                                     par_x03=0.01,
                                     par_x04=-5.00,
                                     par_x05=0.01,
                                     par_x06=0.01,
                                     par_x07=0.01,
                                     par_x08=0.01)
Uparams_defaults = MOHYSE_OST.params(par_x01=20.0,
                                     par_x02=1.0,
                                     par_x03=20.0,
                                     par_x04=5.0,
                                     par_x05=0.5,
                                     par_x06=1.0,
                                     par_x07=1.0,
                                     par_x08=1.0)

hrus_defaults = MOHYSE_OST.hrus(par_x09=0.9039, par_x10=5.6167)
Uhrus_defaults = MOHYSE_OST.hrus(par_x09=15.0,
                                 par_x10=15.0)
Lhrus_defaults = MOHYSE_OST.hrus(par_x09=0.01,
                                 par_x10=0.01)

params = LiteralInput('params', 'Comma separated list of model parameters',
                      abstract='Parameters: ' + ', '.join(params_defaults._fields),
                      data_type='string',
                      default=', '.join(str(p) for p in list(params_defaults)),
                      min_occurs=0)

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

algorithm = LiteralInput('algorithm', 'OSTRICH Algorithm to use to calibrate model parameters',
                         abstract='Optimization algorithm to implement for this calibration run',
                         data_type='string',
                         default='DDS',
                         allowed_values=('DDS', 'SCEUA'),
                         min_occurs=0)

MaxEvals = LiteralInput('MaxEvals', 'Maximum number of model evaluations for the calibration run (budget)',
                        abstract='Maximum number of times OSTRICH can call the hydrological model during the '
                                 'model parameter calibrationn',
                        data_type='integer',
                        default=50,
                        allowed_values=list(range(25001)),
                        min_occurs=0)

CalibrationResults = ComplexOutput('CalibrationResults',
                                   'ObjectiveFunction and calibrated parameters computed by Ostrich',
                                   abstract="Objective Function value after calibration using user-selected "
                                            "function, as well as the calibrated parameter set",
                                   supported_formats=[FORMATS.TEXT],
                                   as_reference=True)


class OstrichMOHYSEProcess(OstrichProcess):
    """
    OSTRICH emulator for the MOHYSE model.

    This process calibrates the MOHYSE model using a OSTRICH emulator. Users need to provide netCDF input
    files for
    rain, snow minimum and maximum temperature as well as potential evapotranspiration. To run diagnostics, observed
    stream flows are also required.
    """
    identifier = 'ostrich-mohyse',
    abstract = 'OSTRICH calibration of RAVEN MOHYSE hydrological model',
    title = '',
    version = '',
    model_cls = MOHYSE_OST,
    tuple_inputs = {'params': MOHYSE_OST.params, 'hrus': MOHYSE_OST.hrus},
    inputs = [wio.ts, algorithm, MaxEvals, params, upperBounds, lowerBounds, wio.start_date, wio.end_date,
              wio.run_name,
              wio.name, wio.area, wio.latitude, wio.longitude, wio.elevation]

    outputs = [CalibrationResults]

    keywords = ["Ostrich", "Calibration", "DDS"],

    def _handler(self, request, response):
        response.update_status('PyWPS process {} started.'.format(self.identifier), 0)
        """
        ts = [e.file for e in request.inputs.pop('ts')]
        algorithm = request.inputs.pop('algorithm')[0].data
        MaxEvals = request.inputs.pop('MaxEvals')[0].data
        params = request.inputs.pop('params')[0].data
        lowerBounds = request.inputs.pop('lowerBounds')[0].data
        upperBounds = request.inputs.pop('upperBounds')[0].data
        start_date = request.inputs.pop('start_date')[0].data
        end_date = request.inputs.pop('end_date')[0].data
        longitude = request.inputs.pop('longitude')[0].data
        run_name = request.inputs.pop('run_name')[0].data
        name = request.inputs.pop('name')[0].data
        area = request.inputs.pop('area')[0].data
        latitude = request.inputs.pop('latitude')[0].data
        elevation = request.inputs.pop('elevation')[0].data

        pdb.set_trace()
        # Write output
        FPath = Path(self.workdir) / 'CalibrationResults.txt'
        with open(FPath, 'w') as f:
            params.to_csv(f, header=['par_x01', 'par_x02', 'par_x03', 'par_x04', 'par_x05',
                                     'par_x06', 'par_x07', 'par_x08', 'par_x09', 'par_x10'])
        response.outputs['CalibrationResults'].file = str(FPath)
        """
