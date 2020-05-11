import six
import datetime as dt
import collections
from pathlib import Path
from xclim.utils import units
from xclim.utils import units2pint

# Can be removed when xclim is pinned above 0.14
units.define("deg_C = degC")

"""
Raven configuration
-------------------

The RV class is used to store Raven parameters for emulated models.

Each model should subclass RV to define the parameters it expects using a namedtuple class. For example::

    class MyModel(RV):
        params = namedtuple('ModelParams', 'x1, x2, x3')
        init = namedtuple('ModelInit', 'i1, i2')
        hru = namedtuple('ModelHRU', 'hru1', hru2')

It can then be instantiated by passing values that will set as default values for each parameter::

    rv = MyModel(params=MyModel.params(1,2,3), init=MyModel.init(0,0), hru=MyModel.hru(4,5), name='basin')

values can then be modified either using attributes or properties::

    rv.name = 'LacVert'
    rv['evaluation_metrics'] = 'LOG_NASH'


Simulation end date and duration are updated automatically when duration, start date or end date are changed.

"""

default_input_variables = ("pr", "rainfall", "prsn", "tasmin", "tasmax", "tas", "evspsbl",
                           "water_volume_transport_in_river_channel")

rain_snow_fraction_options = ("RAINSNOW_DATA", "RAINSNOW_DINGMAN", "RAINSNOW_UBC", "RAINSNOW_HBV", "RAINSNOW_HARDER",
                              "RAINSNOW_HSPF")

evaporation_options = ("PET_CONSTANT", "PET_PENMAN_MONTEITH", "PET_PENMAN_COMBINATION", "PET_PRIESTLEY_TAYLOR",
                       "PET_HARGREAVES", "PET_HARGREAVES_1985", "PET_FROMMONTHLY", "PET_DATA", "PET_HAMON_1961",
                       "PET_TURC_1961", "PET_MAKKINK_1957", "PET_MONTHLY_FACTOR", "PET_MOHYSE", "PET_OUDIN")


class RVFile:

    def __init__(self, fn):
        """Read the content."""
        fn = Path(fn)

        self.stem = fn.with_suffix('').with_suffix('').stem
        self.suffixes = ''.join(fn.suffixes)

        self.ext = ""
        self._store_ext(fn)

        # Whether extension indicates an Ostrich template file.
        self.is_tpl = fn.suffix in ['.tpl', '.txt']

        self.content = ""
        self.content = fn.read_text()

    def _store_ext(self, fn):
        try:
            self.ext = fn.suffixes[0][1:]
        except IndexError as e:
            msg = "\nFile {} does not look like a valid Raven/Ostrich config file.".format(fn)
            raise ValueError(msg) from e

    def rename(self, name):
        self.stem = name

    def write(self, path, **kwds):
        fn = (path / self.stem).with_suffix(self.suffixes)

        content = self.content
        if kwds:
            content = content.format(**kwds)

        fn.write_text(content)
        return fn

    @property
    def tags(self):
        """Return a list of tags within the templates."""
        import re
        pattern = re.compile(r"{([\.\w]+)}")

        return pattern.findall(self.content)


class RV(collections.Mapping):
    """Generic configuration class.

    RV provides two mechanisms to set values, a dictionary-like interface and an object-like interface::

        rv = RV(a=None)
        rv['a'] = 1
        rv.a = 2

    The dictionary like interface only allows the modification of values for existing items, while the object interface
    allows the creation of new attributes::

      rv['c'] = 1

    will raise an AttributeError, while::

      rv.c = 1

    will create a new `c` attribute and assign it the value 1.

    """

    def __init__(self, **kwargs):
        # Set initial default values
        for key, val in kwargs.items():
            setattr(self, key, val)

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        if not hasattr(self, key):
            raise AttributeError('Trying to assign unrecognized object: {}'.format(key))

        setattr(self, key, value)

    def __len__(self):
        return len(self.__dict__)

    def __iter__(self):
        return iter(self.keys())

    def keys(self):
        return (key[1:] if key.startswith('_') else key for key in self.__dict__)

    def items(self):
        for attribute in self.keys():
            yield attribute, getattr(self, attribute)

    def update(self, items, force=False):
        """Update values from dictionary items.

        Parameters
        ----------
        items : dict
          Dictionary of values.
        force : bool
          If True, un-initialized keys can be set.
        """
        if force:
            for key, val in items.items():
                setattr(self, key, val)
        else:
            for key, val in items.items():
                self[key] = val


class RavenNcData(RV):
    pat = """
          :{kind} {raven_name} {site} {runits}
              :ReadFromNetCDF
                 :FileNameNC      {path}
                 :VarNameNC       {var_name}
                 :DimNamesNC      {dimensions}
                 :StationIdx      {index}
                 {time_shift}
                 {linear_transform}
              :EndReadFromNetCDF
          :End{kind}
          """

    _var_names = {'tasmin': "TEMP_MIN",
                  'tasmax': "TEMP_MAX",
                  'tas': "TEMP_AVE",
                  'rainfall': "RAINFALL",
                  'pr': "PRECIP",
                  'prsn': "SNOWFALL",
                  'evspsbl': "PET",
                  'water_volume_transport_in_river_channel': "HYDROGRAPH"
                  }

    _var_runits = {'tasmin': 'deg_C',
                   'tasmax': 'deg_C',
                   'tas': 'deg_C',
                   'pr': "mm/d",
                   'rainfall': "mm/d",
                   'prsn': "mm/d",
                   'evspsbl': "mm/d",
                   'water_volume_transport_in_river_channel': "m3/s"
                   }

    def __init__(self, **kwargs):
        self.var = None
        self.path = None
        self.var_name = None
        self.units = None
        self.scale_factor = None
        self.add_offset = None

        self._time_shift = None
        self._dimensions = None
        self._index = None
        self._linear_transform = None
        self._site = None
        self._kind = None
        self._runits = None
        self._raven_name = None
        self._site = None

        super().__init__(**kwargs)

    @property
    def raven_name(self):
        return self._var_names[self.var]

    @property
    def runits(self):
        return self._var_runits[self.var]

    @property
    def kind(self):
        if self.var == 'water_volume_transport_in_river_channel':
            return "ObservationData"
        else:
            return "Data"

    @property
    def site(self):
        if self.var == 'water_volume_transport_in_river_channel':
            return 1
        else:
            return ""

    @property
    def index(self):
        if self._index is not None:
            return str(self._index + 1)

    @index.setter
    def index(self, value):
        self._index = value

    @property
    def dimensions(self):
        """Return dimensions as Raven expects it:
        - time
        - station time
        - lon lat time
        """
        dims = list(self._dimensions)

        # Move the time dimension at the end
        dims.remove('time')
        dims.append('time')

        return ' '.join(dims)

    @dimensions.setter
    def dimensions(self, value):
        if 'time' not in value:
            raise ValueError("Raven expects a time dimension.")

        self._dimensions = value

        # If there is no spatial dimension, set the index to 0.
        if self.index is None and len(value) == 1:
            self._index = 0

    @property
    def time_shift(self):
        """The fraction describing the time shift, for example to convert UTC to local time."""
        if self._time_shift is not None:
            return f":TimeShift {self._time_shift}"

    @time_shift.setter
    def time_shift(self, value):
        self._time_shift = value

    @property
    def linear_transform(self):
        """A sequence of two values: multiplicative factor and additive offset."""
        lt = self._linear_transform
        if lt is not None or self.scale_factor is not None or self.add_offset is not None:
            slope, intercept = lt or (1, 0)
            sf = 1 if self.scale_factor is None else self.scale_factor
            offset = 0 if self.add_offset is None else self.add_offset
            return ":LinearTransform {:g} {:g}".format(slope * sf, offset * slope + intercept)

    @linear_transform.setter
    def linear_transform(self, value):
        """Set the scale factor and offset."""
        if len(value) == 1:
            value = (value, 0)
        elif len(value) == 2:
            pass
        else:
            raise ValueError("Linear transform takes at most two values: "
                             "a scaling factor and an offset.")

        self._linear_transform = value

    def _check_units(self):
        import warnings
        if units2pint(self.units) != units2pint(self.runits):
            if self._linear_transform is None:
                warnings.warn(f"Units are not what Raven expects for {self.var}.\n"
                              f"Actual: {self.units}\n"
                              f"Expected: {self.runits}\n"
                              f"Make sure to set linear_transform to perform conversion."
                              )

    def __str__(self):
        if self.var is None:
            return ""
        else:
            kwds = {k: v if v is not None else "" for (k, v) in self.items()}
            return self.pat.format(**kwds)


class MonthlyAverage(RV):
    pat = ":MonthlyAve{var}, {data}"

    def __init__(self, var=None, data=None):
        self.var = var
        self.data = data

    def __str__(self):
        if self.var is None:
            return ""
        out = self.pat.format(var=self.var, data=', '.join([str(d) for d in self.data]))
        return out


class RVT(RV):
    def __init__(self, **kwargs):
        self._nc_index = None
        super(RVT, self).__init__(**kwargs)

    @property
    def nc_index(self):
        return self._nc_index

    @nc_index.setter
    def nc_index(self, value):
        for key, val in self.items():
            if isinstance(val, RavenNcData):
                setattr(val, 'index', value)

    def update(self, items, force=False):
        """Update values from dictionary items.

        Parameters
        ----------
        items : dict
          Dictionary of values.
        force : bool
          If True, un-initialized keys can be set.
        """
        if force:
            raise ValueError("Cannot add a new variable at run-time.")
        else:
            for key, val in items.items():
                if isinstance(val, dict):
                    self[key].update(val, force=True)


class RVI(RV):
    def __init__(self, **kwargs):
        self.name = None
        self.area = None
        self.elevation = None
        self.latitude = None
        self.longitude = None
        self.run_index = 0
        self.raven_version = '2.9 rev#177'

        self._run_name = 'run'
        self._start_date = None
        self._end_date = None
        self._now = None
        self._rain_snow_fraction = "RAINSNOW_DATA"
        self._evaporation = None
        self._ow_evaporation = None
        self._duration = 1
        self._time_step = 1.0
        self._evaluation_metrics = 'NASH_SUTCLIFFE RMSE'
        self._suppress_output = False

        super(RVI, self).__init__(**kwargs)

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
        else:
            raise ValueError("Must be datetime")

        if x != dt.datetime(1, 1, 1):
            self._update_duration()

    @property
    def end_date(self):
        return self._end_date

    @end_date.setter
    def end_date(self, x):
        if isinstance(x, dt.datetime):
            self._end_date = x
        else:
            raise ValueError("Must be datetime")

        if x != dt.datetime(1, 1, 1):
            self._update_duration()

    @property
    def duration(self):
        return self._duration

    @duration.setter
    def duration(self, x):
        if isinstance(x, int):
            if x > 0:
                self._duration = x
        else:
            raise ValueError("Must be int")

        if x > 0:
            self._update_end_date()

    @property
    def time_step(self):
        return self._time_step

    @time_step.setter
    def time_step(self, x):
        self._time_step = x

    @property
    def evaluation_metrics(self):
        return self._evaluation_metrics

    @evaluation_metrics.setter
    def evaluation_metrics(self, x):
        if not isinstance(x, six.string_types):
            raise ValueError("Evaluation metrics must be string.")

        for metric in x.split():
            if metric not in {'NASH_SUTCLIFFE', 'LOG_NASH', 'RMSE', 'PCT_BIAS', 'ABSERR', 'ABSMAX', 'PDIFF', 'TMVOL',
                              'RCOEFF', 'NSC', 'KLING_GUPTA'}:
                raise ValueError("{} is not a metric recognized by Raven.")

        self._evaluation_metrics = x

    def _update_duration(self):
        if self.end_date is not None and self.start_date is not None:
            self._duration = (self.end_date - self.start_date).days

    def _update_end_date(self):
        if self.start_date is not None and self.duration is not None:
            self._end_date = self.start_date + dt.timedelta(days=self.duration)

    @property
    def now(self):
        return dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    @property
    def suppress_output(self):
        tag = ":SuppressOutput\n:DontWriteWatershedStorage"
        return tag if self._suppress_output else ""

    @suppress_output.setter
    def suppress_output(self, value):
        if not isinstance(value, bool):
            raise ValueError
        self._suppress_output = value

    @property
    def rain_snow_fraction(self):
        """Rain snow partitioning.
        """
        return self._rain_snow_fraction

    @rain_snow_fraction.setter
    def rain_snow_fraction(self, value):
        """Can be one of

        - RAINSNOW_DATA
        - RAINSNOW_DINGMAN
        - RAINSNOW_UBC
        - RAINSNOW_HBV
        - RAINSNOW_HARDER
        - RAINSNOW_HSPF
        """
        v = value.upper()

        if v in rain_snow_fraction_options:
            self._rain_snow_fraction = v
        else:
            raise ValueError(f"Value should be one of {rain_snow_fraction_options}.")

    @property
    def evaporation(self):
        """Evaporation scheme"""
        return self._evaporation

    @evaporation.setter
    def evaporation(self, value):
        v = value.upper()
        if v in evaporation_options:
            self._evaporation = v
        else:
            raise ValueError(f"Value {v} should be one of {evaporation_options}.")

    @property
    def ow_evaporation(self):
        """Open-water evaporation scheme"""
        return self._ow_evaporation

    @ow_evaporation.setter
    def ow_evaporation(self, value):
        v = value.upper()
        if v in evaporation_options:
            self._ow_evaporation = v
        else:
            raise ValueError(f"Value {v} should be one of {evaporation_options}.")


class RVC(RV):
    def __init__(self):

        # HRU State Variables
        self.surface_water = 0
        self.atmosphere = 0
        self.atmos_precip = 0
        self.ponded_water = 0
        self.soil0 = 0
        self.soil1 = 0
        self.soil2 = 0
        self.soil3 = 0
        self.snow_temp = 0
        self.snow = 0
        self.snow_cover = 0
        self.aet = 0
        self.convolution0 = 0
        self.convolution1 = 0
        self.conv_stor0 = 0
        self.conv_stor1 = 0
        self.conv_stor2 = 0
        self.conv_stor3 = 0
        self.conv_stor4 = 0
        self.conv_stor5 = 0
        self.conv_stor6 = 0
        self.conv_stor7 = 0
        self.conv_stor8 = 0
        self.conv_stor9 = 0
        self.conv_stor10 = 0
        self.conv_stor11 = 0
        self.conv_stor12 = 0
        self.conv_stor13 = 0
        self.conv_stor14 = 0
        self.conv_stor15 = 0
        self.conv_stor16 = 0
        self.conv_stor17 = 0
        self.conv_stor18 = 0
        self.conv_stor19 = 0
        self.conv_stor20 = 0
        self.conv_stor21 = 0
        self.conv_stor22 = 0
        self.conv_stor23 = 0
        self.conv_stor24 = 0
        self.conv_stor25 = 0
        self.conv_stor26 = 0
        self.conv_stor27 = 0
        self.conv_stor28 = 0
        self.conv_stor29 = 0
        self.conv_stor30 = 0
        self.conv_stor31 = 0
        self.conv_stor32 = 0
        self.conv_stor33 = 0
        self.conv_stor34 = 0
        self.conv_stor35 = 0
        self.conv_stor36 = 0
        self.conv_stor37 = 0
        self.conv_stor38 = 0
        self.conv_stor39 = 0
        self.conv_stor40 = 0
        self.conv_stor41 = 0
        self.conv_stor42 = 0
        self.conv_stor43 = 0
        self.conv_stor44 = 0
        self.conv_stor45 = 0
        self.conv_stor46 = 0
        self.conv_stor47 = 0
        self.conv_stor48 = 0
        self.conv_stor49 = 0
        self.conv_stor50 = 0
        self.conv_stor51 = 0
        self.conv_stor52 = 0
        self.conv_stor53 = 0
        self.conv_stor54 = 0
        self.conv_stor55 = 0
        self.conv_stor56 = 0
        self.conv_stor57 = 0
        self.conv_stor58 = 0
        self.conv_stor59 = 0
        self.conv_stor60 = 0
        self.conv_stor61 = 0
        self.conv_stor62 = 0
        self.conv_stor63 = 0
        self.conv_stor64 = 0
        self.conv_stor65 = 0
        self.conv_stor66 = 0
        self.conv_stor67 = 0
        self.conv_stor68 = 0
        self.conv_stor69 = 0
        self.conv_stor70 = 0
        self.conv_stor71 = 0
        self.conv_stor72 = 0
        self.conv_stor73 = 0
        self.conv_stor74 = 0
        self.conv_stor75 = 0
        self.conv_stor76 = 0
        self.conv_stor77 = 0
        self.conv_stor78 = 0
        self.conv_stor79 = 0
        self.conv_stor80 = 0
        self.conv_stor81 = 0
        self.conv_stor82 = 0
        self.conv_stor83 = 0
        self.conv_stor84 = 0
        self.conv_stor85 = 0
        self.conv_stor86 = 0
        self.conv_stor87 = 0
        self.conv_stor88 = 0
        self.conv_stor89 = 0
        self.conv_stor90 = 0
        self.conv_stor91 = 0
        self.conv_stor92 = 0
        self.conv_stor93 = 0
        self.conv_stor94 = 0
        self.conv_stor95 = 0
        self.conv_stor96 = 0
        self.conv_stor97 = 0
        self.conv_stor98 = 0
        self.conv_stor99 = 0

        # Basin state variables
        self.channelstorage = 0
        self.rivuletstorage = 0
        self.qout = [0, 0]
        self.qlat = [0, 0, 0, 0]
        self.qin = 21 * [0, ]

    def initialize(self, path):
        """Set initial conditions based on *solution* output file.

        Parameters
        ----------
        path : str, Path
          Path to `solution.rvc` file.
        """
        with open(path) as f:
            rvc = parse_solution(f.read())

        rvc['HRUStateVariableTable']

    @property
    def hru_state(self):
        pass
    # To be completed

    @property
    def basin_state(self):
        """Return basin state variables."""
        pat = """
              :BasinIndex 1,None
                :ChannelStorage, {self.channelstorage}
                :RivuletStorage, {self.rivuletstorage}
                :Qout,1,{self.qout}
                :Qlat,3,{self.qlat}
                :Qin ,{self.qin}
            """
        return

class Ost(RV):
    def __init__(self, **kwargs):
        self._max_iterations = None
        self._random_seed = None

        super(Ost, self).__init__(**kwargs)

    @property
    def max_iterations(self):
        return self._max_iterations

    @max_iterations.setter
    def max_iterations(self, x):
        if x < 1:
            raise ValueError("Max iteration should be a positive integer: {}".format(x))
        else:
            self._max_iterations = x

    @property
    def random_seed(self):
        if self._random_seed is not None:
            return "RandomSeed {}".format(self._random_seed)
        return ""

    @random_seed.setter
    def random_seed(self, value):
        if value >= 0:
            self._random_seed = value
        else:
            self._random_seed = None


def isinstance_namedtuple(x):
    a = isinstance(x, tuple)
    b = getattr(x, '_fields', None) is not None
    return a and b


def guess_linear_transform(actual, expected):
    """Return RVT compatible dictionary for variable unit transformations.

    Parameters
    ----------
    actual : dict
      The units of each variable.
    expected : dict
      The units expected by Raven.

    Returns
    -------
    dict
      Dictionary keyed by <variable_name>_linear_transform, storing "<scale> <offset>"
      strings used by Raven to transform units.

    """
    # TODO : For precip we also need the frequency to sum over one day.


def parse_solution(rvc):
    """Parse solution file and return dictionary of parameters that can then be used to reinitialize the model.

    Parameters
    ----------
    rvc : str
      Content of a solution.rvc file.
    """
    # Create a generator that will consume lines one by one.
    # tags = ['{i.' + a.lower().replace('[', '').replace(']', '') + '}' for a in atts]
    lines = iter(rvc.splitlines())
    return _parser(lines)


def _parser(lines, indent=""):
    import re
    import itertools
    header_pat = re.compile(r"(\s*):(\w+),?\s*(.*)")

    out = {}
    old_key = None
    for line in lines:
        header = header_pat.match(line)
        if header:
            new_indent, key, value = header.groups()
            if new_indent > indent:
                out[old_key] = _parser(itertools.chain([line,], lines), new_indent)
            elif new_indent < indent:
                return out
            else:
                if key == 'BasinIndex':
                    i, name = value.split(',')
                    out[key] = {i: _parser(lines, new_indent + '  '),
                                'name': name}
                else:
                    out[key] = value.split(',') if ',' in value else value

            old_key = key
        else:
            data = line.split(',')
            i = int(data.pop(0))
            out[i] = list(map(float, data))

    return out
