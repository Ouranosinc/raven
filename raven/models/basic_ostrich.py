"""
Raven model definition

Classes
-------

Ostrich : A generic class that knows how calibrate a model if ostIn.txt and all required files are given.

GR4JCemaneige: The Ostrich calibration of a Raven emulator for GR4J-Cemaneige. Uses template configuration files
whose value can be
automatically filled in.

"""
import raven
from pathlib import Path
from collections import OrderedDict, namedtuple
import os, stat
import subprocess
import tempfile
import csv
import datetime as dt
import six
import xarray as xr
from raven.models import Raven
from raven.models.rv import RV, RVI, isinstance_namedtuple
import numpy as np


def make_executable(fn):
    """Make file executable."""
    st = os.stat(fn)
    os.chmod(fn, st.st_mode | stat.S_IEXEC)


# TODO: Configure this according to the model_path and output_path.
save_best = """#!/bin/bash

set -e

cp ./model/*.rv?  ../best/
cp ./model/output/* ../best/

exit 0
"""

# TODO: Configure this according to raven_cmd, name and output_path.
ostrich_runs_raven = """
#!/bin/bash

set -e

cp ./*.rv? model/

./model/raven ./model/raven-gr4j-salmon -o ./model/output/

exit 0
"""


class Ostrich(Raven):
    """Wrapper for OSTRICH calibration of RAVEN hydrological model

    This class is used to calibrate RAVEN model using OSTRICH from user-provided configuration files. It can also be
    subclassed with
    configuration templates for emulated models, allowing direct calls to the models.

    Usage
    -----
    >>> r = Ostrich('/tmp/testdir')
    >>> r.configure()

    Attributes
    ----------
    conf
      The rv configuration files + Ostrict ostIn.txt
    tpl
      The Ostrich templates

    """
    identifier = 'generic-ostrich'
    txt = RV()

    @property
    def cmd(self):
        """OSTRICH executable path."""
        return self.ostrich_cmd

    @property
    def cmd_path(self):
        """This is the main executable."""
        return self.exec_path

    @property
    def ostrich_cmd(self):
        """OSTRICH executable path."""
        return self.exec_path / 'ostrich'

    @property
    def best_path(self):
        """Path to the best output."""
        return self.exec_path / 'best'

    @property
    def output_path(self):
        """Path to the model outputs and logs."""
        return self.best_path

    def write_save_best(self):
        fn = self.exec_path / 'save_best.sh'
        fn.write_text(save_best)
        make_executable(fn)

    def write_ostrich_runs_raven(self):
        fn = self.exec_path / 'ostrich-runs-raven.sh'
        fn.write_text(ostrich_runs_raven)
        make_executable(fn)

    def setup_model(self, ts, overwrite=False):
        """Create directory structure to store model input files, executable and output results.

        Model configuration files and time series inputs are stored directly in the working directory.

        workdir/  # Created by PyWPS.
           *.rv?
           *.tpl
           ostIn.txt
           model/
           output/

        At each Ostrich loop, configuration files (original and created from templates are copied into model/.

        """
        Raven.setup_model(self, ts, overwrite)

        os.makedirs(str(self.best_path), exist_ok=True)

        self.write_ostrich_runs_raven()
        self.write_save_best()

        # Create symbolic link to executable
        os.symlink(self.ostrich_exec, str(self.cmd))


class GR4JCN_OST(Ostrich):
    """templates = (tuple(Path('ostrich-gr4j-cemaneige').glob("*.tpl")) +
                 tuple(Path('ostrich-gr4j-cemaneige').glob("*.txt")) +
                 tuple(Path('ostrich-gr4j-cemaneige/model').glob("*.rv*")))
    """
    class RVP(RV):
        params = namedtuple('GR4JParams', ('GR4J_X1', 'GR4J_X2', 'GR4J_X3', 'GR4J_X4', 'CEMANEIGE_X1', 'CEMANEIGE_X2'))

    rvp = RVP(params=RVP.params(None, None, None, None, None, None))
    rvt = RV(pr=None, prsn=None, tasmin=None, tasmax=None, evspsbl=None, water_volume_transport_in_river_channel=None)
    rvi = RVI()
    rvh = RV(name=None, area=None, elevation=None, latitude=None, longitude=None)
    rvd = RV(one_minus_CEMANEIGE_X2=None, GR4J_X1_hlf=None)

    def derived_parameters(self):
        self.rvd.GR4J_X1_hlf = self.rvp.params.GR4J_X1 * 1000. / 2.
        self.rvd.one_minus_CEMANEIGE_X2 = 1.0 - self.rvp.params.CEMANEIGE_X2

# class MOHYSE(Raven):
#     templates = tuple((Path(__file__).parent / 'raven-mohyse').glob("*.rv?"))

#     class RVP(RV):
#         params = namedtuple('MOHYSEParams', ', '.join(['par_x{:02}'.format(i) for i in range(1, 9)]))

#     class RVH(RV):
#         hrus = namedtuple('MOHYSEHRU', ('par_x09', 'par_x10'))

#     rvp = RVP(params=RVP.params(*((None, ) * 8)))
#     rvh = RVH(name=None, area=None, elevation=None, latitude=None, longitude=None, hrus=RVH.hrus(None, None))
#     rvt = RV(pr=None, prsn=None, tasmin=None, tasmax=None, evspsbl=None, water_volume_transport_in_river_channel=None)
#     rvi = RVI()
#     rvd = RV(par_rezi_x10=None)

#     def derived_parameters(self):
#         self.rvd['par_rezi_x10'] = 1.0 / self.rvh.hrus.par_x10


# class HMETS(GR4JCN):
#     templates = tuple((Path(__file__).parent / 'raven-hmets').glob("*.rv?"))

#     class RVP(RV):
#         params = namedtuple('HMETSParams', ('GAMMA_SHAPE', 'GAMMA_SCALE', 'GAMMA_SHAPE2', 'GAMMA_SCALE2',
#                                             'MIN_MELT_FACTOR', 'MAX_MELT_FACTOR', 'DD_MELT_TEMP', 'DD_AGGRADATION',
#                                             'SNOW_SWI_MIN', 'SNOW_SWI_MAX', 'SWI_REDUCT_COEFF', 'DD_REFREEZE_TEMP',
#                                             'REFREEZE_FACTOR', 'REFREEZE_EXP', 'PET_CORRECTION',
#                                             'HMETS_RUNOFF_COEFF', 'PERC_COEFF', 'BASEFLOW_COEFF_1',
#                                             'BASEFLOW_COEFF_2', 'TOPSOIL', 'PHREATIC'))

#     rvp = RVP(params=RVP.params(*((None,) * len(RVP.params._fields))))
#     rvt = RV(pr=None, prsn=None, tasmin=None, tasmax=None, evspsbl=None, water_volume_transport_in_river_channel=None)
#     rvi = RVI()
#     rvd = RV(TOPSOIL_m=None, PHREATIC_m=None, SUM_MELT_FACTOR=None, SUM_SNOW_SWI=None, TOPSOIL_hlf=None,
#              PHREATIC_hlf=None)

#     def derived_parameters(self):
#         self.rvd['TOPSOIL_hlf'] = self.rvp.params.TOPSOIL * 0.5
#         self.rvd['PHREATIC_hlf'] = self.rvp.params.PHREATIC * 0.5
#         self.rvd['TOPSOIL_m'] = self.rvp.params.TOPSOIL / 1000.
#         self.rvd['PHREATIC_m'] = self.rvp.params.PHREATIC / 1000.
#         self.rvd['SUM_MELT_FACTOR'] = self.rvp.params.MIN_MELT_FACTOR + self.rvp.params.MAX_MELT_FACTOR
#         self.rvd['SUM_SNOW_SWI'] = self.rvp.params.SNOW_SWI_MIN + self.rvp.params.SNOW_SWI_MAX


# class HBVEC(GR4JCN):
#     templates = tuple((Path(__file__).parent / 'raven-hbv-ec').glob("*.rv?"))

#     class RVP(RV):
#         params = namedtuple('HBVECParams', ('par_x{:02}'.format(i) for i in range(1, 22)))

#     rvp = RVP(params=RVP.params(*((None,) * len(RVP.params._fields))))

#     class RVD(RV):
#         mae = namedtuple('MeanAverageEvap', ('mae_{:02}'.format(i) for i in range(1, 13)))
#         mat = namedtuple('MeanAverageTemp', ('mat_{:02}'.format(i) for i in range(1, 13)))

#     rvd = RVD(one_plus_par_x15=None, par_x11_half=None)

#     rvt = RV(pr=None, prsn=None, tasmin=None, tasmax=None, evspsbl=None,
#              water_volume_transport_in_river_channel=None)

#     rvh = RV(name=None, area=None, elevation=None, latitude=None, longitude=None)

#     def derived_parameters(self):
#         import xarray as xr

#         self.rvd['one_plus_par_x15'] = self.rvp.params.par_x15 + 1.0
#         self.rvd['par_x11_half'] = self.rvp.params.par_x11 / 2.0

#         tasmax = xr.open_dataset(self.rvt.tasmax)[self.rvd.tasmax_var]
#         tasmin = xr.open_dataset(self.rvt.tasmin)[self.rvd.tasmin_var]
#         evap = xr.open_dataset(self.rvt.evspsbl)[self.rvd.evspsbl_var]

#         tas = (tasmax + tasmin) / 2.
#         self.rvd.mat = self.RVD.mat(*tas.groupby('time.month').mean().values)
#         self.rvd.mae = self.RVD.mae(*evap.groupby('time.month').mean().values)
