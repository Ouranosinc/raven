from raven import config
from raven.processes import RavenProcess
from pywps import LiteralInput
from raven.models import MOHYSE
from . import wpsio as wio

# Defaults for this process
params_defaults = MOHYSE.params(par_x01=1.0000,
                                par_x02=0.0468,
                                par_x03=4.2952,
                                par_x04=2.6580,
                                par_x05=0.4038,
                                par_x06=0.0621,
                                par_x07=0.0273,
                                par_x08=0.0453)

hrus_defaults = MOHYSE.hrus(par_x09=0.9039, par_x10=5.6167)

params = LiteralInput('params', 'Comma separated list of model parameters',
                      abstract='Parameters: ' + ', '.join(params_defaults._fields),
                      data_type='string',
                      default=', '.join(str(p) for p in list(params_defaults)),
                      min_occurs=0,
                      max_occurs=config.max_parallel_processes)

hrus = LiteralInput('hrus', 'Comma separated list of HRU parameters',
                    abstract='Parameters: ' + ', '.join(hrus_defaults._fields),
                    data_type='string',
                    default=', '.join(str(p) for p in list(hrus_defaults)),
                    min_occurs=0,
                    max_occurs=config.max_parallel_processes)


class RavenMOHYSEProcess(RavenProcess):
    identifier = 'raven-mohyse'
    abstract = 'MOHYSE hydrological model'
    title = 'TODO'
    version = ''
    model_cls = MOHYSE
    tuple_inputs = {'params': MOHYSE.params, 'hrus': MOHYSE.hrus}

    inputs = [wio.ts, wio.nc_spec, params, hrus, wio.start_date, wio.end_date, wio.duration, wio.run_name,
              wio.name, wio.area, wio.latitude, wio.longitude, wio.elevation, wio.rain_snow_fraction]
