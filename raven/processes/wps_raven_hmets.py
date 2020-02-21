from raven import config
from raven.processes import RavenProcess
from pywps import LiteralInput
from raven.models import HMETS
from . import wpsio as wio

# Defaults for this process
params_defaults = HMETS.params(GAMMA_SHAPE=9.5019,
                               GAMMA_SCALE=0.2774,
                               GAMMA_SHAPE2=6.3942,
                               GAMMA_SCALE2=0.6884,
                               MIN_MELT_FACTOR=1.2875,
                               MAX_MELT_FACTOR=5.4134,
                               DD_MELT_TEMP=2.3641,
                               DD_AGGRADATION=0.0973,
                               SNOW_SWI_MIN=0.0464,
                               SNOW_SWI_MAX=0.1998,
                               SWI_REDUCT_COEFF=0.0222,
                               DD_REFREEZE_TEMP=-1.0919,
                               REFREEZE_FACTOR=2.6851,
                               REFREEZE_EXP=0.3740,
                               PET_CORRECTION=1.0000,
                               HMETS_RUNOFF_COEFF=0.4739,
                               PERC_COEFF=0.0114,
                               BASEFLOW_COEFF_1=0.0243,
                               BASEFLOW_COEFF_2=0.0069,
                               TOPSOIL=310.7211,
                               PHREATIC=916.1947)


params = LiteralInput('params', 'Comma separated list of model parameters',
                      abstract='Parameters: ' + ', '.join(params_defaults._fields),
                      data_type='string',
                      default=', '.join(str(p) for p in list(params_defaults)),
                      min_occurs=0,
                      max_occurs=config.max_parallel_processes)


class RavenHMETSProcess(RavenProcess):
    identifier = 'raven-hmets'
    abstract = 'HMETS hydrological model'
    title = ''
    version = ''
    model_cls = HMETS
    tuple_inputs = {'params': HMETS.params}

    inputs = [wio.ts, wio.nc_spec, params, wio.start_date, wio.end_date, wio.duration, wio.run_name,
              wio.name, wio.area, wio.latitude, wio.longitude, wio.elevation, wio.rain_snow_fraction]
