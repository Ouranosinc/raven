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
from collections import OrderedDict, namedtuple
import os
import subprocess
import csv
import datetime as dt
import six
import xarray as xr
from .rv import RV, RVI

raven_exec = Path(raven.__file__).parent.parent / 'bin' / 'raven'


class Raven:
    """Generic class for the Raven model."""
    identifier = 'generic-raven'
    templates = ()
    rvi = rvp = rvc = rvt = rvh = RV()

    # Output files default names. The actual output file names will be composed of the run_name and the default name.
    _output_fn = {'hydrograph': 'Hydrographs.nc',
                  'storage': 'WatershedStorage.nc',
                  'solution': 'solution.rvc',
                  'diagnostics': 'Diagnostics.csv'}

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

    def __init__(self, workdir):
        self.workdir = Path(workdir)
        self.outputs = {}
        self._name = None

        self._defaults = {}
        self._rvext = ('rvi', 'rvp', 'rvc', 'rvh', 'rvt')

        # The configuration file content is stored in conf.
        self._conf = dict.fromkeys(self._rvext, "")

        # Model parameters - dictionary representation of rv attributes.
        self._parameters = dict.fromkeys(self._rvext, OrderedDict())

        if self.templates:
            self.configure(self.templates)

    @property
    def cmd(self):
        return self.model_path / 'raven'

    @property
    def model_path(self):
        return self.workdir / 'model'

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, x):
        if self._name is None:
            self._name = x
        elif x != self._name:
            raise UserWarning("Model configuration name changed.")

    @property
    def output_path(self):
        return self.workdir / 'output'

    @property
    def parameters(self):
        return {ext: OrderedDict(getattr(self, ext).to_dict()) for ext in self._rvext}

    @property
    def rv(self):
        """Return a dictionary of the configuration files."""
        return {ext: self._conf[ext] for ext in self._rvext}

    @property
    def rvobjs(self):
        """Return a generator for (ext, rv object)."""
        return {ext: getattr(self, ext) for ext in self._rvext}

    def configure(self, fns):
        """Read configuration files."""
        for fn in fns:
            name, ext = self.split_ext(fn)
            self.name = name

            if ext not in self._rvext:
                raise ValueError('rv contains unrecognized configuration file keys : {}.'.format(ext))

            self._conf[ext] = open(fn).read()

    def assign(self, key, value):
        """Assign parameter to rv object that has a key with the same name."""
        for ext, obj in self.rvobjs.items():
            try:
                obj[key] = value
            except AttributeError:
                pass

    def derived_parameters(self):
        return {}

    def _dump_rv(self, ts):
        """Write configuration files to disk."""

        for ext, txt in self.rv.items():
            fn = self.model_path / (self.name + '.' + ext)

            with open(fn, 'w') as f:

                # Write parameters into template. Don't write by default as some configuration files might have
                # meaningless {}.
                if self.parameters[ext]:
                    params = self.parameters[ext]

                    # Include the variable names to the parameters
                    if ext == 'rvt':
                        files, var_names = self._assign_files(ts, self.rvt.keys())
                        params.update(files)
                        params.update(var_names)

                    txt = txt.format(**params)

                f.write(txt)

    def setup_model(self, ts, overwrite=False):
        """Create directory structure to store model input files, executable and output results.

        Model configuration files and time series inputs are stored directly in the working directory.

        workdir/  # Created by PyWPS. Is considered the model path.
           model/
           output/

        """
        import shutil

        # Create subdirectory
        os.makedirs(self.output_path, exist_ok=overwrite)
        os.makedirs(self.model_path, exist_ok=overwrite)

        # Compute derived parameters
        self.derived_parameters()

        # Write configuration files in model directory
        self._dump_rv(ts)

        # Create symbolic link to input files
        for fn in ts:
            os.symlink(fn, self.model_path / Path(fn).name)

        # Create symbolic link to executable
        os.symlink(raven_exec, self.cmd)

    def run(self, ts, overwrite=False, **kwds):
        """Run the model.

        Parameters
        ----------
        ts : sequence
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
        >>> r = Raven('/tmp/test')
        >>> r.configure(rvi=<path to template>, rvp=...}
        >>> r.run(rvp={'param1': val1, ...}, rcv={...})
        """
        # Update parameter objects
        for key, val in kwds.items():
            obj = getattr(self, key)
            if isinstance(val, dict):
                obj.update(val)
            else:
                obj.values = val

        if self.rvi:
            self.handle_date_defaults(ts)

        self.setup_model(map(Path, ts), overwrite)

        # Run the model
        subprocess.call(map(str, [self.cmd, self.model_path / self.name, '-o', self.output_path]))

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
        try:
            ds = xr.open_mfdataset(fns)
        except OSError:
            raise NotImplementedError
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
    def hydrograph(self):
        import xarray as xr
        return xr.open_dataset(self.outputs['hydrograph'])

    @property
    def storage(self):
        import xarray as xr
        return xr.open_dataset(self.outputs['storage'])

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

    @staticmethod
    def split_ext(fn):
        """Return the name and rv key of the configuration file."""
        if isinstance(fn, six.string_types):
            fn = Path(fn)

        return (fn.stem, fn.suffix[1:])


class GR4JCemaneige(Raven):
    templates = tuple((Path(__file__).parent / 'raven-gr4j-cemaneige').glob("*.rv?"))

    class RVP(RV):
        params = namedtuple('GR4JParams', ('GR4J_X1', 'GR4J_X2', 'GR4J_X3', 'GR4J_X4', 'CEMANEIGE_X1', 'CEMANEIGE_X2'))

    rvp = RVP(params=RVP.params(None, None, None, None, None, None), one_minus_CEMANEIGE_X2=None)
    rvc = RV(GR4J_X1_hlf=None)
    rvt = RV(pr=None, prsn=None, tasmin=None, tasmax=None, evspsbl=None, water_volume_transport_in_river_channel=None)
    rvi = RVI()
    rvh = RV(name=None, area=None, elevation=None, latitude=None, longitude=None)

    def derived_parameters(self):
        self.rvc.GR4J_X1_hlf = self.rvp.params.GR4J_X1 * 1000. / 2.
        self.rvp.one_minus_CEMANEIGE_X2 = 1.0 - self.rvp.params.CEMANEIGE_X2


class MOHYSE(Raven):
    templates = tuple((Path(__file__).parent / 'raven-mohyse').glob("*.rv?"))

    class RVP(RV):
        params = namedtuple('MOHYSEParams', ', '.join(['par_x{:02}'.format(i) for i in range(1, 9)]))

    class RVH(RV):
        hrus = namedtuple('MOHYSEHRU', ('par_x09', 'par_x10'))

    rvp = RVP(params=RVP.params(*((None, ) * 8)))
    rvh = RVH(name=None, area=None, elevation=None, latitude=None, longitude=None, hrus=RVH.hrus(None, None),
              par_rezi_x10=None)
    rvt = RV(pr=None, prsn=None, tasmin=None, tasmax=None, evspsbl=None, water_volume_transport_in_river_channel=None)
    rvi = RVI()

    def derived_parameters(self):
        self.rvh['par_rezi_x10'] = 1.0 / self.rvh.hrus.par_x10


class HMETS(GR4JCemaneige):
    templates = tuple((Path(__file__).parent / 'raven-hmets').glob("*.rv?"))

    class RVP(RV):
        params = namedtuple('HMETSParams', ('GAMMA_SHAPE', 'GAMMA_SCALE', 'GAMMA_SHAPE2', 'GAMMA_SCALE2',
                                            'MIN_MELT_FACTOR', 'MAX_MELT_FACTOR', 'DD_MELT_TEMP', 'DD_AGGRADATION',
                                            'SNOW_SWI_MIN', 'SNOW_SWI_MAX', 'SWI_REDUCT_COEFF', 'DD_REFREEZE_TEMP',
                                            'REFREEZE_FACTOR', 'REFREEZE_EXP', 'PET_CORRECTION',
                                            'HMETS_RUNOFF_COEFF', 'PERC_COEFF', 'BASEFLOW_COEFF_1',
                                            'BASEFLOW_COEFF_2', 'TOPSOIL', 'PHREATIC'))

    rvp = RVP(params=RVP.params(*((None,) * len(RVP.params._fields))), TOPSOIL_m=None, PHREATIC_m=None,
              SUM_MELT_FACTOR=None, SUM_SNOW_SWI=None)
    rvc = RV(TOPSOIL_hlf=None, PHREATIC_hlf=None)
    rvt = RV(pr=None, prsn=None, tasmin=None, tasmax=None, evspsbl=None, water_volume_transport_in_river_channel=None)
    rvi = RVI()

    def derived_parameters(self):
        self.rvc['TOPSOIL_hlf'] = self.rvp.params.TOPSOIL * 0.5
        self.rvc['PHREATIC_hlf'] = self.rvp.params.PHREATIC * 0.5
        self.rvp['TOPSOIL_m'] = self.rvp.params.TOPSOIL / 1000.
        self.rvp['PHREATIC_m'] = self.rvp.params.PHREATIC / 1000.
        self.rvp['SUM_MELT_FACTOR'] = self.rvp.params.MIN_MELT_FACTOR + self.rvp.params.MAX_MELT_FACTOR
        self.rvp['SUM_SNOW_SWI'] = self.rvp.params.SNOW_SWI_MIN + self.rvp.params.SNOW_SWI_MAX
