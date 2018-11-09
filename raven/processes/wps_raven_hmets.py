from collections import OrderedDict as Odict
from raven.processes import RavenProcess
from pywps import LiteralInput
from raven.models import HMETS
from . import wpsio as wio

# Defaults for this process
param_defaults = Odict([('GAMMA_SHAPE', 9.5019),
                   ('GAMMA_SCALE', 0.2774),
                   ('GAMMA_SHAPE2', 6.3942),
                   ('GAMMA_SCALE2', 0.6884),
                   ('MIN_MELT_FACTOR', 1.2875),
                   ('MAX_MELT_FACTOR', 6.7009),
                   ('DD_MELT_TEMP', 2.3641),
                   ('DD_AGGRADATION', .0973),
                   ('SNOW_SWI_MIN', 0.0464),
                   ('SNOW_SWI_MAX', 0.2462),
                   ('SWI_REDUCT_COEFF', 0.0222),
                   ('DD_REFREEZE_TEMP', -1.0919),
                   ('REFREEZE_FACTOR', 2.6851),
                   ('REFREEZE_EXP', 0.3740),
                   ('PET_CORRECTION', 1.0),
                   ('HMETS_RUNOFF_COEFF', 0.4739),
                   ('PERC_COEFF', 0.0114),
                   ('BASEFLOW_COEFF_1', 0.0243),
                   ('BASEFLOW_COEFF_2', 0.0069),
                   ('TOPSOIL', 0.3107),
                   ('PHREATIC', 0.9162)])

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


class RavenHMETSProcess(RavenProcess):
    identifier = 'raven-hmets'
    abstract = 'HMETS hydrological model'
    title = ''
    version = ''
    model_cls = HMETS
    param_arrays = ['params', 'init']

    inputs = [wio.ts, params, wio.start_date, wio.end_date, wio.duration, init, wio.run_name,
              wio.name, wio.area, wio.latitude, wio.longitude, wio.elevation]


