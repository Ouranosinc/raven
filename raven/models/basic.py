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
import subprocess
import csv
import six
import datetime as dt

raven_exec = Path(raven.__file__).parent.parent / 'bin' / 'raven'

class Raven:
    """Generic class for the Raven model."""
    name = 'generic-raven'

    # Output files default names. The actual output file names will be composed of the run_name and the default name.
    _output_fn = {'hydrograph': 'Hydrographs.nc',
                  'storage': 'WatershedStorage.nc',
                  'solution': 'solution.rvc',
                  'diagnostics': 'Diagnostics.csv'}

    def __init__(self, workdir):
        self.workdir = Path(workdir)

        self._defaults = {}
        self._outputs = {}
        self._rvext = ('rvi', 'rvp', 'rvc', 'rvh', 'rvt')

        # Configuration files dictionary (can also be template)
        self._configuration = dict.fromkeys(self._rvext, '')

        # Model parameters
        self.parameters = dict.fromkeys(self._rvext, OrderedDict())

    @property
    def model_path(self):
        return self.workdir / 'model'

    @property
    def output_path(self):
        return self.workdir / 'output'

    @property
    def cmd(self):
        return self.model_path / 'raven'

    def configure(self, **kwds):
        for key, fn in kwds.items():
            if key not in self._rvext:
                raise ValueError('rv contains unrecognized configuration file keys.')

            self._configuration[key] = Path(fn)
            setattr(self, key, open(fn).read())

        self.name = fn.stem

    @property
    def configuration(self):
        """Return rv file paths."""
        return self._configuration

    def _dump_rv(self):
        """Write configuration files to disk."""
        for key, val in self.configuration.items():
            fn = self.model_path / val.name
            with open(fn, 'w') as f:
                f.write(getattr(self, key).format(**self.parameters[key]))

    def setup_model(self, fns, overwrite=False):
        """Create directory structure to store model input files, executable and output results.

        Model configuration files and time series inputs are stored directly in the working directory.

        workdir/  # Created by PyWPS. Is considered the model path.
           model/
           output/

        """


        # Create subdirectory
        os.makedirs(self.output_path, exist_ok=overwrite)
        os.makedirs(self.model_path, exist_ok=overwrite)

        # Write configuration files in model directory
        self._dump_rv()

        # Create symbolic link to input files
        for fn in fns:
            os.symlink(fn, self.model_path / fn.name)

        # Create symbolic link to executable
        os.symlink(raven_exec, self.cmd)

    def run(self, fns, overwrite=False, **kwds):
        """Run the model.

        Parameters
        ----------
        fns : sequence
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
        for key, val in kwds.items():
            self.parameters[key].update(val)

        self.setup_model(fns, overwrite)

        # Run the model
        subprocess.call(map(str, [self.cmd, self.model_path / self.name, '-o', self.output_path]))

        # Store output file names in dict
        for key in self._output_fn.keys():
            self._outputs[key] = self._get_output(key)

    __call__ = run

    def _get_output(self, key):
        """Match actual output files to known expected files.

        Return a dictionary of file paths for each expected input.
        """
        fn = self._output_fn[key]
        files = self.output_path.glob('*' + fn)
        if len(files) == 0:
            raise IOError("No output files for {}".format(fn))
        if len(files) > 1:
            raise IOError("Multiple matching files found for {}.".format(fn))

        return files[0].absolute()

    @property
    def hydrograph(self):
        import xarray as xr
        return xr.open_dataset(self._outputs['hydrograph'])

    @property
    def storage(self):
        import xarray as xr
        return xr.open_dataset(self._outputs['storage'])

    @property
    def diagnostics(self):
        with open(self._outputs['diagnostics']) as f:
            reader = csv.reader(f.readlines())
            header = reader.next()
            content = reader.next()

            out = dict(zip(header, content))
            out.pop('')
            return out


class RVI:

    def __init__(self, **kwds):
        for key, val in kwds:
            if val is not None:
                setattr(self, key, val)
            else:
                setattr(self, '_'+key, None)

    @property
    def run_name(self):
        return self._run_name

    @run_name.setter
    def run_name(self, x):
        if isinstance(x, six.string_types):
            self._run_name = x
        else:
            raise ValueError("Must be string")

    @property
    def start_date(self):
        return self._start_date

    @start_date.setter
    def start_date(self, x):
        if isinstance(x, dt.datetime):
            self._start_date = x
            self._update_duration()
        else:
            raise ValueError("Must be datetime")

    @property
    def end_date(self):
        return self._end_date

    @end_date.setter
    def end_date(self, x):
        if isinstance(x, dt.datetime):
            self._end_date = x
            self._update_duration()
        else:
            raise ValueError("Must be datetime")

    @property
    def duration(self):
        if isinstance(int):
            return self._duration
        else:
            raise ValueError("Must be int")

    @duration.setter
    def duration(self, x):
        self._duration = x
        self._update_end_date()

    @property
    def time_step(self):
        return self._time_step

    @time_step.setter
    def time_step(self, x):
        self._time_step = x

    def _update_duration(self):
        if self.end_date is not None and self.start_date is not None:
            self.duration = (self.end_date - self.start_date).days

    def _update_end_date(self):
        if self.start_date is not None and self.duration is not None:
            self.end_date = self.start_date + dt.timedelta(days=1)



