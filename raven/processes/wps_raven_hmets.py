from collections import OrderedDict as Odict

from raven.processes import RavenGR4JCemaNeigeProcess

# from pywps import LiteralInput

# Model defaults
defaults = Odict(
    rvi=dict(run_name=None, Start_Date=None, End_Date=None, Duration=None, TimeStep=1.0,
             EvaluationMetrics='NASH_SUTCLIFFE RMSE'),
    rvp=Odict([('GAMMA_SHAPE', None),
               ('GAMMA_SCALE', None),
               ('GAMMA_SHAPE2', None),
               ('GAMMA_SCALE2', None),
               ('MIN_MELT_FACTOR', None),
               ('MAX_MELT_FACTOR', None),
               ('DD_MELT_TEMP', None),
               ('DD_AGGRADATION', None),
               ('SNOW_SWI_MIN', None),
               ('SNOW_SWI_MAX', None),
               ('SNOW_SWI', 0.5),
               ('SWI_REDUCT_COEFF', None),
               ('DD_REFREEZE_TEMP', None),
               ('REFREEZE_FACTOR', None),
               ('REFREEZE_EXP', None),
               ('PET_CORRECTION', None),
               ('HMETS_RUNOFF_COEFF', None),
               ('PERC_COEFF', None),
               ('BASEFLOW_COEFF_1', None),
               ('BASEFLOW_COEFF_2', None),
               ('TOPSOIL', None),
               ('PHREATIC', None)]),
    rvc=Odict(SOIL_0=None, SOIL_1=None),
    rvh=dict(NAME=None, AREA=None, ELEVATION=None, LATITUDE=None, LONGITUDE=None),
    rvt=dict(pr=None, prsn=None, tasmin=None, tasmax=None, evspsbl=None, water_volume_transport_in_river_channel=None)
)


# Defaults for this process


class RavenHMETSProcess(RavenGR4JCemaNeigeProcess):
    identifier = 'raven-hmets'
    abstract = 'HMETS hydrological model'
    title = ''
    defaults = defaults
    pdefaults = Odict([('GAMMA_SHAPE', 9.5019),
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
