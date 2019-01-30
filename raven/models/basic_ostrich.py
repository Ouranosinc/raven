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


save_best = """#!/bin/bash

set -e

echo "saving input files for the best solution found ..."

if [ ! -e best ] ; then
    mkdir best
fi

cp ./model/*.rv?  ../best/
cp ./model/output/* ../best/

exit 0
"""

ostric_runs_raven = """
#!/bin/bash

set -e

echo "Running raven ..."

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
    ostrich_exec = raven.ostrich_exec

    @property
    def ostrich_cmd(self):
        """OSTRICH executable path."""
        return self.run_path / 'ostrich'

    @property
    def best_path(self):
        """Path to the best output."""
        return self.run_path / 'best'

    @property
    def output_path(self):
        """Path to the model outputs and logs."""
        return self.best_path

    @property
    def tpl(self):
        """Dictionary of Ostrich template files content keyed by extension."""
        return self._tpl

    @property
    def rv(self):
        """Dictionary of the configuration files."""
        return {ext: self._conf.get(ext, self._tpl.get(ext, "")) for ext in self._rvext}

    def configure(self, fns):
        """Read configuration files."""

        for fn in fns:
            print(">>>>>>>>>>> fn = ", fn)
            name, ext = self.split_ext(fn)

            if ext == 'tpl':  # Ostrich template
                name, ext = self.split_ext(name)
                store = self._tpl
            else:
                store = self._conf

            self._name = name

            if ext not in self._rvext:
                raise ValueError('rv contains unrecognized configuration file keys : {}.'.format(ext))

            store[ext] = open(str(fn)).read()

    def assign(self, key, value):
        """
        Assign parameter to rv object that has a key with the same name.
        For example, {run_name}, {longitude}, {elevation}.
        """
        assigned = False
        for ext, obj in self.rvobjs.items():
            if hasattr(obj, key):
                att = getattr(obj, key)
                # If the object is a namedtuple, we get its class and try to instantiate it with the values passed.
                if isinstance_namedtuple(att) and isinstance(value, (list, tuple, np.ndarray)):
                    p = getattr(getattr(self, ext.upper()), key)(*value)
                    setattr(obj, key, p)
                else:
                    setattr(obj, key, value)
                assigned = True

        if not assigned:
            raise AttributeError("No configuration key named {}".format(key))

    def derived_parameters(self):
        """Subclassed by emulators. Defines model parameters that are a function of other parameters."""
        return

    def _dump_rv(self, ts):
        """Write configuration files to disk."""

        # Merge all parameter dictionaries
        params = {}
        for key, val in self.configuration.items():
            params.update(val)

        for ext in self._rvext:
            if ext in self.conf:
                txt = self.conf[ext]

                if ext == 'txt':
                    fn = str(self.run_path / 'ostIn.txt')
                else:
                    fn = str(self.model_path / (self.name + '.' + ext))

            elif ext in self.tpl:
                txt = self.tpl[ext]
                fn = str(self.run_path / (self.name + '.' + ext + '.tpl'))

            else:
                continue

            with open(fn, 'w') as f:
                # Write parameters into template.
                if params:
                    txt = txt.format(**params)

                f.write(txt)

    def write_save_best(self):
        fn = self.run_path / 'save_best.sh'
        with open(fn, 'w') as f:
            f.write(save_best)
        make_executable(fn)

    def write_ostrich_runs_raven(self):
        fn = self.run_path / 'ostrich-runs-raven.sh'
        with open(fn, 'w') as f:
            f.write(ostric_runs_raven)
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
        import shutil

        if self.output_path.exists():
            if overwrite:
                shutil.rmtree(str(self.workdir))
            else:
                raise IOError(
                    "Directory already exists. Either set overwrite to `True` or create a new model instance.")

        # Create subdirectory
        os.makedirs(str(self.run_path))
        os.makedirs(str(self.model_path), exist_ok=True)
        os.makedirs(str(self.output_path), exist_ok=True)
        os.makedirs(str(self.best_path), exist_ok=True)

        self.write_ostrich_runs_raven()
        self.write_save_best()

        # Match the input files
        files, var_names = self._assign_files(ts, self.rvt.keys())

        self.rvt.update(files, force=True)
        self.rvd.update(var_names, force=True)

        # Compute derived parameters
        self.derived_parameters()

        # needs to fill in {run_name}, {start_date}, {duration}, {time_step}, etc. --> model/model/*.rv*  OR
        # model/*.rv*.tpl
        # needs to fill in information {calib_alg}, {calib_budget}, {calib_ranges} --> model/ostIn.txt

        # Write configuration files in model directory
        self._dump_rv(ts)  # RAVEN specific

        # Create symbolic link to input files
        for fn in ts:
            os.symlink(str(fn), str(self.model_path / Path(fn).name))

        # Raven executable
        os.symlink(self.raven_exec, str(self.raven_cmd))

        # Create symbolic link to executable
        os.symlink(self.ostrich_exec, str(self.ostrich_cmd))

    def run(self, ts, overwrite=False, **kwds):
        """Run the model.

        Parameters
        ----------
        ts : path or sequence
          Sequence of input file paths. Symbolic links to those files will be created in the model directory.
        overwrite : bool
          Whether or not to overwrite existing model and output files.
        **kwds : dict
          Raven parameters used to fill configuration file templates.

        Create a work directory with a model/ and output/ subdirectories, write the configuration files in model/ and
        launch the Raven executable. If the configuration files are templates, values can be formatted by passing
        dictionaries keyed by their extension.

        Example
        -------
        >>> r = Ostrich()
        >>> r.configure(rvi=<path to template>, rvp=...}
        >>> r.run(ts, start_date=dt.datetime(2000, 1, 1), area=1000, X1=67)
        """
        import shutil

        if isinstance(ts, (six.string_types, Path)):
            ts = [ts, ]

        # Update parameter objects
        for key, val in kwds.items():

            if key in self._rvext:
                obj = getattr(self, key)
                if isinstance(val, dict):
                    obj.update(val)
                elif isinstance(val, RV):
                    setattr(self, key, val)
                else:
                    raise ValueError("A dictionary or an RV instance is expected to update the values "
                                     "for {}.".format(key))
            else:

                self.assign(key, val)

            print("KJDHSSSSSSSSSSS: ", key, "   ", val)

        if self.rvi:
            self.handle_date_defaults(ts)

        self.setup_model(tuple(map(Path, ts)), overwrite)

        # Run the model (OSTRICH needs to be executed in directory of executable)
        prevdir = os.getcwd()
        os.chdir(self.run_path)
        subprocess.call(map(str, [self.ostrich_cmd]))
        os.chdir(prevdir)

        # Store output file names in dict
        for key in self._output_fn.keys():
            self.outputs[key] = self._get_output(key)

    __call__ = run

    def _assign_files(self, fns, variables):
        """Find for each variable the file storing it's data and the name of the netCDF variable.

        Parameters
        ----------
        fns : sequence
          Paths to netCDF files.
        variables : sequence
          Names of the variables to look for. Specify their CF standard name, a dictionary of
          alternative names will be used for the lookup.

        Returns
        -------
        files : dict
          A dictionary keyed by variable storing the file storing each variable.
        variables : dict
          A dictionary keyed by variable_var storing the variable name within the netCDF file.
        """
        files = {}
        var_names = {}

        for fn in fns:
            if '.nc' in fn.suffix:
                with xr.open_dataset(fn) as ds:
                    for var in variables:
                        for alt_name in self._variable_names[var]:
                            if alt_name in ds.data_vars:
                                files[var] = fn
                                var_names[var + '_var'] = alt_name
                                break

        for var in variables:
            if var not in files.keys():
                raise ValueError("{} not found in files.".format(var))

        return files, var_names

    def _get_output(self, key):
        """Match actual output files to known expected files.

        Return a dictionary of file paths for each expected input.
        """
        fn = self._output_fn[key]
        files = list(self.output_path.glob('*' + fn))

        if len(files) == 0:
            raise UserWarning("No output files for {}".format(fn))

        if len(files) > 1:
            raise IOError("Multiple matching files found for {}.".format(fn))

        return str(files[0].absolute())

    @staticmethod
    def start_end_date(fns):
        """Return the common starting and ending date and time of netCDF files.

        Parameters
        ----------
        fns : sequence
          Sequence of netCDF file names for forcing data.

        Returns
        -------
        start : datetime
          The first datetime of the forcing files.
        end : datetime
          The last datetime of the forcing files.
        """

        ds = xr.open_mfdataset(fns)
        return ds.indexes['time'][0], ds.indexes['time'][-1]

    def handle_date_defaults(self, ts):

        # Get start and end date from file
        start, end = self.start_end_date(ts)

        rvi = self.rvi
        if rvi.start_date == dt.datetime(1, 1, 1):
            rvi.start_date = start

        else:
            if rvi.end_date == dt.datetime(1, 1, 1):
                rvi.end_date = end

    @property
    def q_sim(self):
        """Return a view of the hydrograph time series.

        This view will be overwritten by successive calls to `run`. To make a copy of this DataArray that will
        persist in memory, use `q_sim.copy(deep=True)`.
        """
        return self.hydrograph.q_sim

    @property
    def hydrograph(self):
        """Return a view of the current output file.

        If the model is run multiple times, hydrograph will point to the latest version. To store the results of
        multiple runs, either create different model instances or explicitly copy the file to another disk location.
        """
        with xr.open_dataset(self.outputs['hydrograph']) as ds:
            return ds

    @property
    def storage(self):
        with xr.open_dataset(self.outputs['storage']) as ds:
            return ds

    @property
    def diagnostics(self):
        with open(self.outputs['diagnostics']) as f:
            reader = csv.reader(f.readlines())
            header = next(reader)
            content = next(reader)

            out = dict(zip(header, content))
            out.pop('')

        for key, val in out.items():
            if 'DIAG' in key:
                out[key] = float(val)

        return out

    @property
    def tags(self):
        """Return a list of tags within the templates."""
        import re
        pattern = re.compile(r"{(\w+)}")

        out = {}
        if self.templates:
            for key, conf in self.rv.items():
                out[key] = pattern.findall(conf)

        return out


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
