"""
Raven model definition

Classes
-------

Raven : A generic class that knows how to launch the model from completed rv files.
Ostrich:

"""
import csv
import datetime as dt
import os
import shutil
import stat
import subprocess
import tempfile
from collections import OrderedDict
from pathlib import Path

import numpy as np
import six
import xarray as xr

import raven
from .rv import RVFile, RV, RVI, isinstance_namedtuple, Ost, RavenNcData, parse_solution


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

    # Dictionary of potential variable names, keyed by CF standard name.
    # http://cfconventions.org/Data/cf-standard-names/60/build/cf-standard-name-table.html
    # PET is the potential evapotranspiration, while evspsbl is the actual evap.
    # TODO: Check we're not mixing precip and rainfall.
    _variable_names = {'tasmin': ['tasmin', 'tmin'],
                       'tasmax': ['tasmax', 'tmax'],
                       'tas': ['tas', 't2m'],
                       'rainfall': ['rainfall', 'rain'],
                       'pr': ['pr', 'precip', 'prec', 'precipitation', 'tp'],
                       'prsn': ['prsn', 'snow', 'snowfall', 'solid_precip'],
                       'evspsbl': ['pet', 'evap', 'evapotranspiration'],
                       'water_volume_transport_in_river_channel': ['qobs', 'discharge', 'streamflow', 'dis']
                       }

    # Expected units (pint-compatible)
    _units = {'tasmin': "degC",
              'tasmax': "degC",
              'tas': "degC",
              'pr': "mm/d",
              'rainfall': "mm/d",
              'prsn': "mm/d",
              'evspsbl': "mm/d",
              'water_volume_transport_in_river_channel': "m**3/s"
              }

    _parallel_parameters = ['params', 'hru_state', 'basin_state', 'nc_index', 'name', 'area', 'elevation', 'latitude',
                            'longitude', 'region_id', 'hrus']

    def __init__(self, workdir=None):
        """Initialize the RAVEN model.

        Parameters
        ----------
        workdir : str, Path
          Directory for the model configuration and outputs. If None, a temporary directory will be created.
        """
        workdir = workdir or tempfile.mkdtemp()
        self._rvs = []

        self.rvi = RV()
        self.rvp = RV()
        self.rvc = RV()
        self.rvt = RV()
        self.rvh = RV()
        self.rvd = RV()  # rvd is for derived parameters

        self.workdir = Path(workdir)
        self.ind_outputs = {}  # Individual files for all simulations
        self.outputs = {}  # Aggregated files
        self.singularity = False  # Set to True to launch Raven with singularity.
        self.raven_exec = raven.raven_exec
        self.raven_simg = raven.raven_simg
        self.ostrich_exec = raven.ostrich_exec
        self._name = None
        self._defaults = {}
        self.rvfiles = {}

        # Configuration file extensions + rvd for derived parameters.
        self._rvext = self._rvext + ('rvd',)

        # For subclasses where the configuration file templates are known in advance.
        if self.templates:
            self.configure(self.templates)

        # Directory logic
        # Top directory inside workdir. This is where Ostrich and its config and templates are stored.
        self.model_dir = 'model'  # Path to the model configuration files.
        self.final_dir = 'final'
        self.output_dir = 'output'

        self.exec_path = self.workdir / 'exec'
        self.final_path = self.workdir / self.final_dir
        self._psim = 0
        self._pdim = None  # Parallel dimension (either params or nbasins)

    @property
    def output_path(self):
        return self.model_path / self.output_dir

    @property
    def model_path(self):
        return self.exec_path / self.model_dir / "p{:02}".format(self.psim)

    @property
    def raven_cmd(self):
        """Path to the Raven executable."""
        return self.model_path / 'raven'

    @property
    def version(self):
        import re
        out = subprocess.check_output([self.raven_exec, ], input=b'\n')
        match = re.search(r"Version (\S+) ", out.decode('utf-8'))
        if match:
            return match.groups()[0]
        else:
            raise AttributeError("Version not found: {}".format(out))

    @property
    def psim(self):
        return self._psim

    @psim.setter
    def psim(self, value):
        if not isinstance(value, int):
            raise ValueError
        if isinstance(self.rvi, RVI):
            self.rvi.run_index = value
        self._psim = value

    @property
    def cmd(self):
        """This is the main executable."""
        return self.raven_cmd

    @property
    def bash_cmd(self):
        """Bash command arguments."""
        return [self.cmd, self.name, '-o', str(self.output_path)]

    @property
    def singularity_cmd(self):
        """Run Singularity container."""
        return ["singularity", "run", "--bind", "{}:/data".format(self.model_path), "--bind",
                "{}:/data_out:rw".format(self.output_path), self.raven_simg, self.name]

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
        self._name = x

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
                if rvf.ext.startswith('rv'):
                    setattr(self, 'name', rvf.stem)
                    self.rvfiles[rvf.ext] = rvf
                elif rvf.ext == 'txt':
                    self.rvfiles[rvf.stem] = rvf
                else:
                    raise ValueError

    def assign(self, key, value):
        """Assign parameter to rv object that has a key with the same name."""
        assigned = False
        for ext, obj in self.rvobjs.items():
            if hasattr(obj, key):
                att = getattr(obj, key)

                # If att is a namedtuple, we get its class and try to instantiate it with the values passed.
                if isinstance_namedtuple(att) and isinstance(value, (list, tuple, np.ndarray)):
                    p = att.__class__(*value)
                    setattr(obj, key, p)
                # If att is a RavenNcData, we expect a dict
                elif isinstance(att, RavenNcData):
                    att.update(value)
                else:
                    setattr(obj, key, value)
                assigned = True

        if not assigned:
            raise AttributeError("No configuration key named {}".format(key))

    def derived_parameters(self):
        """Subclassed by emulators. Defines model parameters that are a function of other parameters."""
        return

    def _dump_rv(self):
        """Write configuration files to disk."""
        params = self.parameters

        for rvf in self.rvfiles.values():
            p = self.exec_path if rvf.is_tpl else self.model_path
            if rvf.stem == 'OstRandomNumbers' and isinstance(self.txt, Ost) and self.txt.random_seed == "":
                continue
            fn = rvf.write(p, **params)
            self._rvs.append(fn)

    def setup(self, overwrite=False):
        """Create directory structure to store model input files, executable and output results.

        Model configuration files and time series inputs are stored directly in the working directory.

        workdir/  # Created by PyWPS. Is considered the model path.
           model/
           output/

        """
        if overwrite:
            if self.model_path.exists():
                shutil.rmtree(str(self.exec_path))
            if self.final_path.exists():
                shutil.rmtree(str(self.final_path))

        # Create general subdirectories
        if not self.exec_path.exists():
            os.makedirs(str(self.exec_path))  # workdir/exec
        if not self.final_path.exists():
            os.makedirs(str(self.final_path))  # workdir/final

    def setup_model_run(self, ts):
        """Create directory structure to store model input files, executable and output results.

        Parameters
        ----------
        ts : sequence
          Paths to input forcing files.
        index : int
          Run index.
        """
        # Create configuration information from input files
        ncvars = self._assign_files(ts)
        self.rvt.update(ncvars)
        self.check_units()
        self.check_inputs()

        # Compute derived parameters
        self.derived_parameters()

        # Write configuration files in model directory
        if not self.model_path.exists():
            os.makedirs(self.model_path)
            os.makedirs(self.output_path)
        self._dump_rv()

        # Create symbolic link to input files
        for fn in ts:
            if not (self.model_path / Path(fn).name).exists():
                os.symlink(str(fn), str(self.model_path / Path(fn).name))

        # Create symbolic link to Raven executable
        if not self.raven_cmd.exists():
            os.symlink(self.raven_exec, str(self.raven_cmd))

        # Shell command to run the model
        if self.singularity:
            cmd = self.singularity_cmd
        else:
            cmd = self.bash_cmd

        return cmd

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
        >>> r.configure(rvi='path to template', rvp='...'}
        >>> r.run(ts, start_date=dt.datetime(2000, 1, 1), area=1000, X1=67)

        """
        if isinstance(ts, (six.string_types, Path)):
            ts = [ts, ]

        # Case for potentially parallel parameters
        pdict = {}
        for p in self._parallel_parameters:
            a = kwds.pop(p, None)

            if p in ['params', 'hrus', 'basin_state', 'hru_state'] and a is not None:
                pdict[p] = np.atleast_2d(a)
            else:
                pdict[p] = np.atleast_1d(a)

        # Number of parallel loops is dictated by the number of parameters or nc_index.
        nloops = max(len(pdict['params']), len(pdict['nc_index']))
        if nloops > 1:
            self._pdim = 'nbasins' if len(pdict['nc_index']) > 1 else 'params'

        for key, val in pdict.items():
            if len(val) not in [1, nloops]:
                raise ValueError("Parameter {} has incompatible dimension: {}. "
                                 "Should be 1 or {}.".format(key, len(val), nloops))

        # Resize parallel parameters to the largest size
        for key, val in pdict.items():
            if len(val) == 1:
                pdict[key] = val.repeat(nloops, axis=0)

        # Update non-parallel parameter objects
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
            self.set_calendar(ts)

        # Loop over parallel parameters
        procs = []
        for self.psim in range(nloops):
            for key, val in pdict.items():
                if val[self.psim] is not None:
                    self.assign(key, val[self.psim])

            cmd = self.setup_model_run(tuple(map(Path, ts)))
            procs.append(subprocess.Popen(cmd, cwd=self.cmd_path, stdout=subprocess.PIPE))

        return procs

    def __call__(self, ts, overwrite=False, **kwds):
        self.setup(overwrite)
        procs = self.run(ts, overwrite, **kwds)

        for proc in procs:
            proc.wait()
            # Julie: For debugging
            # for line in iter(proc.stdout.readline, b''):
            #    print(line)
        try:
            self.parse_results()

        except UserWarning as e:
            err = self.parse_errors()
            msg = """
        **************************************************************
        Path : {dir}
        **************************************************************
        {err}
        """.format(dir=self.cmd_path, err=err)
            print(msg)
            raise e

    def resume(self, solution=None):
        """Set the initial state to the state at the end of the last run.

        Parameters
        ----------
        solution : str, Path
          Path to solution file. If None, will use solution from last model run if any.

        Note that the starting date of the next run should be identical to the end date of the previous run.
        """
        if solution is None:
            fn = self.outputs['solution']
        else:
            fn = solution

        rvc = RVFile(fn)
        rvc.rename(self.name)
        self.rvfiles['rvc'] = rvc

    def parse_results(self, path=None):
        """Store output files in the self.outputs dictionary."""
        # Output files default names. The actual output file names will be composed of the run_name and the default
        # name.
        path = path or self.exec_path

        patterns = {'hydrograph': '*Hydrographs.nc',
                    'storage': '*WatershedStorage.nc',
                    'solution': '*solution.rvc',
                    'diagnostics': '*Diagnostics.csv'
                    }

        for key, pattern in patterns.items():
            # There are no diagnostics if a streamflow time series is not provided.
            try:
                fns = self._get_output(pattern, path=path)
            except UserWarning as exc:
                if key != 'diagnostics':
                    raise exc

            fns.sort()
            self.ind_outputs[key] = fns
            self.outputs[key] = self._merge_output(fns, pattern[1:])

        self.outputs['rv_config'] = self._merge_output(self.rvs, 'rv.zip')

    def _merge_output(self, files, name):
        """Merge multiple output files into one if possible, otherwise return a list of files.
        """
        import zipfile

        # If there is only one file, return its name directly.
        if len(files) == 1:
            return files[0]

        # Otherwise try to create a new file aggregating all files.
        outfn = self.final_path / name

        if name.endswith('.nc') and not isinstance(self, raven.models.RavenMultiModel):
            ds = [xr.open_dataset(fn) for fn in files]
            try:
                # We aggregate along the pdim dimensions.
                out = xr.concat(ds, self._pdim, data_vars='different')
                out.to_netcdf(outfn)
                return outfn
            except (ValueError, KeyError):
                pass

        # Let's zip the files that could not be merged.
        outfn = outfn.with_suffix('.zip')

        # Find the lower file parts level at which there are differences among files.
        i = get_diff_level(files)

        # Try to create a zip file
        with zipfile.ZipFile(outfn, 'w') as f:
            for fn in files:
                len(fn.parts)
                f.write(fn, arcname=fn.relative_to(Path(*fn.parts[:i])))

        return outfn

    def parse_errors(self):
        files = self._get_output('Raven_errors.txt', self.exec_path)
        out = ''
        for f in files:
            out += f.read_text()
        return out

    def _assign_files(self, fns):
        """Find for each variable the file storing it's data and the name of the netCDF variable.

        Parameters
        ----------
        fns : sequence
          Paths to netCDF files.

        Returns
        -------
        dict
          A dictionary keyed by variable storing the `RavenNcData` instance storing each variable's configuration
          information.
        """
        ncvars = {}
        for fn in fns:
            if '.nc' in fn.suffix:
                with xr.open_dataset(fn) as ds:
                    for var, alt_names in self._variable_names.items():
                        # Check that the emulator is expecting that variable.
                        if var not in self.rvt.keys():
                            continue

                        # Check if any alternate variable name is in the file.
                        for alt_name in alt_names:
                            if alt_name in ds.data_vars:
                                ncvars[var] = dict(var=var,
                                                   path=fn,
                                                   var_name=alt_name,
                                                   dimensions=ds[alt_name].dims,
                                                   units=ds[alt_name].attrs.get("units"),
                                                   )

                                break
        return ncvars

    def _get_output(self, pattern, path):
        """Match actual output files to known expected files.

        Return a dictionary of file paths for each expected input.
        """
        files = list(path.rglob(pattern))

        # TODO: Fix this. Raven won't have rvi.suppress_output is initialized with existing configuration files.
        if len(files) == 0 and not self.rvi.suppress_output:
            raise UserWarning("No output files for {} in {}.".format(pattern, path))

        return [f.absolute() for f in files]

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

        ds = xr.open_mfdataset(fns, combine="by_coords")
        return ds.indexes['time'][0], ds.indexes['time'][-1]

    @staticmethod
    def get_calendar(fns):
        """Return the calendar."""
        ds = xr.open_mfdataset(fns, combine="by_coords")
        return ds.time.encoding["calendar"]

    def set_calendar(self, ts):
        """Set the calendar in the RVI configuration."""
        self.rvi.calendar = self.get_calendar(ts)

    def handle_date_defaults(self, ts):
        # Get start and end date from file
        start, end = self.start_end_date(ts)

        rvi = self.rvi
        if rvi.start_date in [None, dt.datetime(1, 1, 1)]:
            rvi.start_date = start

        else:
            if rvi.end_date in [None, dt.datetime(1, 1, 1)]:
                rvi.end_date = end

    @property
    def rvs(self):
        return self._rvs

    @property
    def q_sim(self):
        """Return a view of the hydrograph time series.

        This view will be overwritten by successive calls to `run`. To make a copy of this DataArray that will
        persist in memory, use `q_sim.copy(deep=True)`.
        """
        if isinstance(self.hydrograph, list):
            return [h.q_sim for h in self.hydrograph]

        return self.hydrograph.q_sim

    @property
    def hydrograph(self):
        """Return a view of the current output file.

        If the model is run multiple times, hydrograph will point to the latest version. To store the results of
        multiple runs, either create different model instances or explicitly copy the file to another disk location.
        """
        if self.outputs['hydrograph'].suffix == '.nc':
            return xr.open_dataset(self.outputs['hydrograph'])
        elif self.outputs['hydrograph'].suffix == '.zip':
            return [xr.open_dataset(fn) for fn in self.ind_outputs['hydrograph']]
        else:
            raise ValueError

    @property
    def storage(self):
        if self.outputs['storage'].suffix == '.nc':
            return xr.open_dataset(self.outputs['storage'])
        elif self.outputs['storage'].suffix == '.zip':
            return [xr.open_dataset(fn) for fn in self.ind_outputs['storage']]
        else:
            raise ValueError

    @property
    def solution(self):
        if self.outputs['solution'].suffix == ".rvc":
            return parse_solution(self.outputs['solution'].read_text())

    @property
    def diagnostics(self):
        diag = []
        for fn in self.ind_outputs['diagnostics']:
            with open(fn) as f:
                reader = csv.reader(f.readlines())
                header = next(reader)
                content = next(reader)

                out = dict(zip(header, content))
                out.pop('')

            for key, val in out.items():
                if 'DIAG' in key:
                    out[key] = float(val)
            diag.append(out)

        return diag if len(diag) > 1 else diag[0]

    @property
    def tags(self):
        """Return a list of tags within the templates."""
        out = []
        for rvf in self.rvfiles.values():
            out.extend(rvf.tags)

        return out

    @staticmethod
    def split_ext(fn):
        """Return the name and rv key of the configuration file."""
        if isinstance(fn, six.string_types):
            fn = Path(fn)

        return fn.stem, fn.suffix[1:]

    def check_units(self):
        """Check that the input file units match expectations."""
        for var, nc in self.rvt.items():
            if isinstance(nc, RavenNcData) and nc.var is not None:
                nc._check_units()

    def check_inputs(self):
        """Check that necessary variables are defined."""
        has_file = set([key for key, val in self.rvt.items() if val is not None])
        vars = list(self.rvt.keys())

        for var in vars:
            if var not in has_file and var != "nc_index":
                if var in ['tasmin', 'tasmax'] and 'tas' in has_file:
                    pass  # This is OK
                if var == 'tas' and has_file.issuperset(['tasmin', 'tasmax']):
                    pass
                elif var in ['prsn']:
                    pass  # Ok, can be guessed from temp ?
                elif var in ['evspsbl']:
                    pass  # Ok, can be computed by Oudin ?
                elif var in ['water_volume_transport_in_river_channel']:
                    pass  # Ok, not strictly necessary for simulations ?
                else:
                    raise ValueError("{} not found in files.".format(var))


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

    @property
    def model_path(self):
        return self.exec_path / self.model_dir

    @staticmethod
    def _allowed_extensions():
        return Raven._allowed_extensions() + ('txt',)

    @property
    def ostrich_cmd(self):
        """OSTRICH executable path."""
        return self.exec_path / 'ostrich'

    @property
    def cmd(self):
        """OSTRICH executable path."""
        return self.ostrich_cmd

    @property
    def cmd_path(self):
        """This is the main executable."""
        return self.exec_path

    @property
    def proc_path(self):
        """Path to Ostrich parallel process directory."""
        return self.exec_path / 'processor_0'  # /'model' / 'output' ?

    def write_save_best(self):
        fn = self.exec_path / 'save_best.sh'
        fn.write_text(save_best)
        make_executable(fn)

    def write_ostrich_runs_raven(self):
        fn = self.exec_path / 'ostrich-runs-raven.sh'
        fn.write_text(ostrich_runs_raven.format(name=self.name))
        make_executable(fn)

    def setup(self, overwrite=False):
        """Create directory structure to store model input files, executable and output results.

        Model configuration files and time series inputs are stored directly in the working directory.

        workdir/  # Created by PyWPS.
           *.rv?
           *.tpl
           ostIn.txt
           model/
           model/output/
           best/

        At each Ostrich loop, configuration files (original and created from templates are copied into model/.

        """
        Raven.setup(self, overwrite)

        os.makedirs(str(self.final_path), exist_ok=True)

        self.write_ostrich_runs_raven()
        self.write_save_best()

        # Create symbolic link to executable
        os.symlink(self.ostrich_exec, str(self.cmd))

    def parse_results(self):
        """Store output files in the self.outputs dictionary."""
        # Output files default names. The actual output file names will be composed of the run_name and the default
        # name.
        Raven.parse_results(self, path=self.final_path)

        patterns = {'params_seq': 'OstModel?.txt',
                    'calibration': 'OstOutput?.txt',
                    }

        # Store output file names in dict
        for key, pattern in patterns.items():
            self.outputs[key] = self._get_output(pattern, path=self.exec_path)[0]

        self.outputs['calibparams'] = ', '.join(map(str, self.calibrated_params))

    def parse_errors(self):
        try:
            raven_err = self._get_output('OstExeOut.txt', path=self.exec_path)[0].read_text()
        except UserWarning:  # Read in processor_0 directory instead.
            try:
                raven_err = self._get_output('OstExeOut.txt', path=self.proc_path)[0].read_text()
            except UserWarning:
                raven_err = ''

        try:
            ost_err = self._get_output('OstErrors?.txt', path=self.exec_path)[0].read_text()
        except UserWarning:  # Read in processor_0 directory instead.
            ost_err = self._get_output('OstErrors?.txt', path=self.proc_path)[0].read_text()

        return "{}\n{}".format(ost_err, raven_err)

    def parse_optimal_parameter_set(self):
        """Return dictionary of optimal parameter set."""
        import re
        txt = open(self.outputs['calibration']).read()
        ops = re.search(r'.*Optimal Parameter Set(.*?)\n{2}', txt, re.DOTALL).groups()[0]

        p = re.findall(r'(\w+)\s*:\s*([\S]+)', ops)
        return OrderedDict((k, float(v)) for k, v in p)

    def ost2raven(self, ops):
        """Return model parameters.

        Note
        ----
        This method should be subclassed by emulators for which Ostrich has different parameters than the original
        Raven model.
        """
        if hasattr(self, 'params'):
            n = len(self.params._fields)
            pattern = 'par_x{}' if n < 8 else 'par_x{:02}'
            names = [pattern.format(i + 1) for i in range(n)]
            return self.params(*(ops[n] for n in names))
        else:
            return ops.values()

    @property
    def calibrated_params(self):
        """The dictionary of optimal parameters estimated by Ostrich."""
        ops = self.parse_optimal_parameter_set()
        return self.ost2raven(ops)

    @property
    def obj_func(self):
        return np.loadtxt(self.outputs['params_seq'], skiprows=1)[-1, 1]

    @property
    def optimized_parameters(self):
        """These are the raw parameters returned by Ostrich."""
        return np.loadtxt(self.outputs['params_seq'], skiprows=1)[-1, 2:]


def get_diff_level(files):
    """Return the lowest hierarchical file parts level at which there are differences among file paths."""

    for i, parts in enumerate(zip(*[f.parts for f in files])):
        if len(set(parts)) > 1:
            return i


def make_executable(fn):
    """Make file executable."""
    st = os.stat(fn)
    os.chmod(fn, st.st_mode | stat.S_IEXEC)


# TODO: Configure this according to the model_path and output_path.
save_best = """#!/bin/bash

set -e

cp ./model/*.rv?  ../../final/
cp ./model/output/* ../../final/

exit 0
"""

# TODO: Configure this according to raven_cmd, name and output_path.
ostrich_runs_raven = """
#!/bin/bash

set -e

cp ./*.rv? model/

./model/raven ./model/{name} -o ./model/output/

exit 0
"""
