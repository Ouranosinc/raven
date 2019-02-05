from .wps_ostrich import OstrichProcess
from raven.models import HMETS_OST
from . import wpsio as wio
import logging
from pywps import FORMATS, LiteralInput, ComplexOutput
import pdb
from pathlib import Path

LOGGER = logging.getLogger("PYWPS")

"""
Notes
-----
The configuration files for a OSTRICH calibration of the HMETS model and in models/ostrich-hmets.
All parameters that could potentially be user-defined are tagged using {}. These tags need to be replaced by
actual values before the model is launched.
"""

params_defaults = HMETS_OST.params(GAMMA_SHAPE=9.5019,
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

Lparams_defaults = HMETS_OST.params(GAMMA_SHAPE=0.3,
                                    GAMMA_SCALE=0.01,
                                    GAMMA_SHAPE2=0.5,
                                    GAMMA_SCALE2=0.15,
                                    MIN_MELT_FACTOR=0.0,
                                    MAX_MELT_FACTOR=0.0,
                                    DD_MELT_TEMP=-2.0,
                                    DD_AGGRADATION=0.01,
                                    SNOW_SWI_MIN=0.0,
                                    SNOW_SWI_MAX=0.01,
                                    SWI_REDUCT_COEFF=0.005,
                                    DD_REFREEZE_TEMP=-5.0,
                                    REFREEZE_FACTOR=0.0,
                                    REFREEZE_EXP=0.0,
                                    PET_CORRECTION=0.0,
                                    HMETS_RUNOFF_COEFF=0.0,
                                    PERC_COEFF=0.00001,
                                    BASEFLOW_COEFF_1=0.0,
                                    BASEFLOW_COEFF_2=0.00001,
                                    TOPSOIL=0.0,
                                    PHREATIC=0.0)

Uparams_defaults = HMETS_OST.params(GAMMA_SHAPE=20.0,
                                    GAMMA_SCALE=5.0,
                                    GAMMA_SHAPE2=13.0,
                                    GAMMA_SCALE2=1.5,
                                    MIN_MELT_FACTOR=20.0,
                                    MAX_MELT_FACTOR=20.0,
                                    DD_MELT_TEMP=3.0,
                                    DD_AGGRADATION=0.2,
                                    SNOW_SWI_MIN=0.1,
                                    SNOW_SWI_MAX=0.3,
                                    SWI_REDUCT_COEFF=0.1,
                                    DD_REFREEZE_TEMP=2.0,
                                    REFREEZE_FACTOR=5.0,
                                    REFREEZE_EXP=1.0,
                                    PET_CORRECTION=3.0,
                                    HMETS_RUNOFF_COEFF=1.0,
                                    PERC_COEFF=0.02,
                                    BASEFLOW_COEFF_1=0.1,
                                    BASEFLOW_COEFF_2=0.01,
                                    TOPSOIL=0.5,
                                    PHREATIC=2.0)

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


class OstrichHMETSProcess(OstrichProcess):
    """
    OSTRICH emulator for the HMETS model.

    This process calibrates the HMETS model using a OSTRICH emulator. Users need to provide netCDF input
    files for
    rain, snow minimum and maximum temperature as well as potential evapotranspiration. To run diagnostics, observed
    stream flows are also required.
    """
    identifier = 'ostrich-hmets'
    abstract = 'OSTRICH calibration of RAVEN HMETS hydrological model'
    title = ''
    version = ''
    model_cls = HMETS_OST
    tuple_inputs = {'lowerBounds': HMETS_OST.params,
                    'upperBounds': HMETS_OST.params}
    inputs = [wio.ts, lowerBounds, upperBounds, wio.algorithm, wio.MaxEvals, wio.start_date, wio.end_date,
              wio.duration, wio.run_name, wio.name, wio.area, wio.latitude, wio.longitude, wio.elevation]

    keywords = ["Ostrich", "Calibration", "DDS"]
