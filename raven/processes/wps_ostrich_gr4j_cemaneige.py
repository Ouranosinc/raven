from .wps_ostrich import OstrichProcess
from raven.models import GR4JCN_OST
from . import wpsio as wio
import logging
from pywps import FORMATS, LiteralInput, ComplexOutput
import pdb
from pathlib import Path

LOGGER = logging.getLogger("PYWPS")

"""
Notes
-----
The configuration files for a OSTRICH calibration of the GR4J-Cemaneige model and in models/ostrich-gr4j-cemaneige.
All parameters that could potentially be user-defined are tagged using {}. These tags need to be replaced by
actual values before the model is launched.
"""

params_defaults = GR4JCN_OST.params(GR4J_X1=0.529,
                                        GR4J_X2=-3.396,
                                        GR4J_X3=407.29,
                                        GR4J_X4=1.072,
                                        CEMANEIGE_X1=16.9,
                                        CEMANEIGE_X2=0.947)

Uparams_defaults = GR4JCN_OST.params(GR4J_X1=0.9,
                                         GR4J_X2=0.,
                                         GR4J_X3=500.,
                                         GR4J_X4=1.1,
                                         CEMANEIGE_X1=20.,
                                         CEMANEIGE_X2=1.)

Lparams_defaults = GR4JCN_OST.params(GR4J_X1=0.1,
                                         GR4J_X2=-5.,
                                         GR4J_X3=100.,
                                         GR4J_X4=1.,
                                         CEMANEIGE_X1=10,
                                         CEMANEIGE_X2=0.1)

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


class OstrichGR4JCemaNeigeProcess(OstrichProcess):
    """
    OSTRICH emulator for the GR4J-Cemaneige model.

    This process calibrates the GR4J-Cemaneige model using a OSTRICH emulator. Users need to provide netCDF input
    files for
    rain, snow minimum and maximum temperature as well as potential evapotranspiration. To run diagnostics, observed
    stream flows are also required.
    """
    identifier = 'ostrich-gr4j-cemaneige',
    abstract = 'OSTRICH calibration of RAVEN GR4J + CEMANEIGE hydrological model',
    title = '',
    version = '',
    model_cls = GR4JCN_OST,
    tuple_inputs = {'params': GR4JCN_OST.params},
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
            params.to_csv(f, header=['X1', 'X2', 'X3', 'X4', 'X5', 'X6'])
        response.outputs['CalibrationResults'].file = str(FPath)
        """
