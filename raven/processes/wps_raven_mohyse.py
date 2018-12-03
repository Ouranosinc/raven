from collections import OrderedDict as Odict
from raven.processes import RavenProcess
from pywps import LiteralInput
from raven.models import MOHYSE
from . import wpsio as wio

# Defaults for this process
param_defaults = Odict([('par_x01', 1.0000),
                        ('par_x02', 0.0468),
                        ('par_x03', 4.2952),
                        ('par_x04', 2.6580),
                        ('par_x05', 0.4038),
                        ('par_x06', 0.0621),
                        ('par_x07', 0.0273),
                        ('par_x08', 0.0453),
                        ('par_x09', 0.9039),
                        ('par_x10', 5.6167)])


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


class RavenMOHYSEProcess(RavenProcess):
    identifier = 'raven-mohyse'
    abstract = 'MOHYSE hydrological model'
    title = ''
    version = ''
    model_cls = MOHYSE
    param_arrays = ['params', 'init']

    inputs = [wio.ts, params, wio.start_date, wio.end_date, wio.duration, init, wio.run_name,
              wio.name, wio.area, wio.latitude, wio.longitude, wio.elevation]
