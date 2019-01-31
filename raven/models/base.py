"""
Raven model definition

Classes
-------

Raven : A generic class that knows how to launch the model from completed rv files.





GR4JCemaneige: The Raven emulator for GR4J-Cemaneige. Uses template configuration files whose value can be
automatically filled in.

"""
import raven
from pathlib import Path
from collections import OrderedDict
import os
import stat
import subprocess
import tempfile
import csv
import datetime as dt
import six
import xarray as xr
from .rv import RVFile, RV, isinstance_namedtuple
import numpy as np


class Raven:
    """RAVEN hydrological model wrapper

    This class is used to run the RAVEN model from user-provided configuration files. It can also be subclassed with
    configuration templates for emulated models, allowing direct calls to the models.

    Usage
    -----
    >>> r = Raven('/tmp/testdir')
    >>> r.configure()

    """
    identifier = 'generic-raven'
    templates = ()

    # Allowed configuration file extensions
    _rvext = ('rvi', 'rvp', 'rvc', 'rvh', 'rvt')

    rvi = rvp = rvc = rvt = rvh = rvd = RV()  # rvd is for derived parameters

    # Output files default names. The actual output file names will be composed of the run_name and the default name.
    _output_fn = {'hydrograph': 'Hydrographs.nc',
                  'storage': 'WatershedStorage.nc',
                  'solution': 'solution.rvc',
                  'diagnostics': 'Diagnostics.csv',
                  }

    # Dictionary of potential variable names, keyed by CF standard name.
    # http://cfconventions.org/Data/cf-standard-names/60/build/cf-standard-name-table.html
    # PET is the potential evapotranspiration, while evspsbl is the actual evap.
    _variable_names = {'tasmin': ['tasmin', 'tmin'],
                       'tasmax': ['tasmax', 'tmax'],
                       'pr': ['pr', 'precip', 'prec', 'rain', 'rainfall', 'precipitation'],
                       'prsn': ['prsn', 'snow', 'snowfall', 'solid_precip'],
                       'evspsbl': ['pet', 'evap', 'evapotranspiration'],
                       'water_volume_transport_in_river_channel': ['qobs', 'discharge', 'streamflow']
                       }

    def __init__(self, workdir=None):
        """Initialize the RAVEN model.

        Parameters
        ----------
        workdir : str, Path
          Directory for the model configuration and outputs. If None, a temporary directory will be created.
        """
        workdir = workdir or tempfile.mkdtemp()
        self.workdir = Path(workdir)
        self.outputs = {}
        self.raven_exec = raven.raven_exec
        self.ostrich_exec = raven.ostrich_exec
        self._name = None
        self._defaults = {}
        self.rvfiles = []

        # Configuration file extensions + rvd for derived parameters.
        self._rvext = self._rvext + ('rvd', )

        # For subclasses where the configuration file templates are known in advance.
        if self.templates:
            self.configure(self.templates)

        # Directory logic

        # Top directory inside workdir. This is where Ostrich and its config and templates are stored.
        self.exec_path = self.workdir / 'exec'

        # Path to the Raven executable and configuration files.
        self.model_path = self.exec_path / 'model'
        self.raven_cmd = self.model_path / 'raven'

    @property
    def output_path(self):
        """Path to the model outputs and logs."""
        return self.model_path / 'output'

    @property
    def version(self):
        import re
        out = subprocess.check_output([self.raven_exec, ])
        match = re.search(r"Version (\S+) ", out.decode('utf-8'))
        if match:
            return match.groups()[0]
        else:
            raise AttributeError("Version not found: {}".format(out))

    @property
    def cmd(self):
        """This is the main executable."""
        return self.raven_cmd

    @property
    def cmd_path(self):
        """This is the main executable."""
        return self.model_path

    @property
    def name(self):
        """Name of the model configuration."""
        return self._name

    @name.setter
    def name(self, x):
        if self._name is None:
            self._name = x
        elif x != self._name:
            raise UserWarning("Model configuration name changed.")
        try:
            if self.rvi.run_name is None:
                self.rvi.run_name = x
        except AttributeError:
            pass

    @property
    def configuration(self):
        """Configuration dictionaries."""
        return {ext: OrderedDict(getattr(self, ext).items()) for ext in self._rvext}

    @property
    def parameters(self):
        """Dictionary storing all parameters."""
        params = {}
        for key, val in self.configuration.items():
            params.update(val)
        return params

    @property
    def rvobjs(self):
        """Generator for (ext, rv object)."""
        return {ext: getattr(self, ext) for ext in self._rvext}

    def configure(self, fns):
        """Read configuration files."""
        for fn in fns:
            rvf = RVFile(fn)
            if rvf.ext not in self._rvext + ('txt',):
                raise ValueError('rv contains unrecognized configuration file keys : {}.'.format(rvf.ext))
            else:
                if rvf.ext != 'txt':
                    setattr(self, 'name', rvf.stem)
                self.rvfiles.append(rvf)

    def assign(self, key, value):
        """Assign parameter to rv object that has a key with the same name."""
        assigned = False
        for ext, obj in self.rvobjs.items():
            if hasattr(obj, key):
                att = getattr(obj, key)

                # If the object is a namedtuple, we get its class and try to instantiate it with the values passed.
                if isinstance_namedtuple(att) and isinstance(value, (list, tuple, np.ndarray)):
                    p = getattr(self, key)(*value)
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
        params = self.parameters

        for rvf in self.rvfiles:
            p = self.exec_path if rvf.is_tpl else self.model_path
            rvf.write(p, **params)

    def setup_model(self, ts, overwrite=False):
        """Create directory structure to store model input files, executable and output results.

        Model configuration files and time series inputs are stored directly in the working directory.

        workdir/  # Created by PyWPS. Is considered the model path.
           model/
           output/

        """
        import shutil

        if self.output_path.exists():
            if overwrite:
                shutil.rmtree(str(self.workdir))
            else:
                raise IOError(
                    "Directory already exists. Either set overwrite to `True` or create a new model instance.")

        # Create subdirectory
        os.makedirs(str(self.exec_path))
        os.makedirs(str(self.model_path), exist_ok=True)
        os.makedirs(str(self.output_path), exist_ok=True)

        # Match the input files
        files, var_names = self._assign_files(ts, self.rvt.keys())
        self.rvt.update(files, force=True)
        self.rvd.update(var_names, force=True)

        # Compute derived parameters
        self.derived_parameters()

        # Write configuration files in model directory
        self._dump_rv(ts)

        # Create symbolic link to input files
        for fn in ts:
            os.symlink(str(fn), str(self.model_path / Path(fn).name))

        # Create symbolic link to Raven executable
        os.symlink(self.raven_exec, str(self.raven_cmd))

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
        >>> r = Raven()
        >>> r.configure(rvi=<path to template>, rvp=...}
        >>> r.run(ts, start_date=dt.datetime(2000, 1, 1), area=1000, X1=67)
        """
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

        if self.rvi:
            self.handle_date_defaults(ts)

        self.setup_model(tuple(map(Path, ts)), overwrite)

        # Run the model
        try:
            cmd = ['./' + self.cmd.stem, self.name, '-o', str(self.output_path)]
            proc = subprocess.Popen(cmd, cwd=self.cmd_path)
            proc.wait()

        except Exception as e:
            msg = (' '.join(map(str, cmd)))
            print("Executed: \n {}\n in `{}`".format(msg, self.cmd_path))
            raise e

        try:
            # Store output file names in dict
            for key, pattern in self._output_fn.items():
                self.outputs[key] = str(self._get_output(pattern))

        except UserWarning as e:
            print("Work directory: ", self.exec_path)
            msg = self._get_error_message()
            print(msg)
            raise e


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

    def _get_error_message(self):
        return self._get_output('Raven_errors.txt').read_text()

    def _get_output(self, pattern):
        """Match actual output files to known expected files.

        Return a dictionary of file paths for each expected input.
        """
        files = list(self.output_path.glob(pattern))

        if len(files) == 0:
            raise UserWarning("No output files for {}.".format(pattern))

        if len(files) > 1:
            raise IOError("Multiple matching files found for {}.".format(pattern))

        return files[0].absolute()

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
        out = []
        for rvf in self.rvfiles:
            out.extend(rvf.tags)

        return out

    @staticmethod
    def split_ext(fn):
        """Return the name and rv key of the configuration file."""
        if isinstance(fn, six.string_types):
            fn = Path(fn)

        return (fn.stem, fn.suffix[1:])


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
    _rvext = ('rvi', 'rvp', 'rvc', 'rvh', 'rvt', 'txt')
    txt = RV()

    @staticmethod
    def _allowed_extensions():
        return Raven._allowed_extensions() + ('txt', )

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

    def _get_error_message(self):
        raven_err = self._get_output('OstExeOut.txt').read_text()
        ost_err = self._get_output('OstErrors?.txt').read_text()
        return "{}\n{}".format(ost_err, raven_err)


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
