import six
import datetime as dt
import collections


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


class RV(collections.abc.Mapping):
    """Generic configuration class."""

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

    def to_dict(self):
        """Return dictionary of content and namedtuple fields."""
        out = {}
        for key, val in self.items():
            if hasattr(val, '_fields'):
                out.update(val._asdict())
            else:
                out[key] = val
        return out


class RVI(RV):
    """Configuration class for rvp, rvh, rvc"""
    def __init__(self, **kwargs):
        self.name = None
        self.area = None
        self.elevation = None
        self.latitude = None
        self.longitude = None

        self._run_name = None
        self._start_date = None
        self._end_date = None
        self._duration = 1
        self._time_step = 1.0
        self._evaluation_metrics = 'NASH_SUTCLIFFE RMSE'

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
            self._update_end_date()

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


def isinstance_namedtuple(x):
    a = isinstance(x, tuple)
    b = isinstance(getattr(x, '__dict__', None), collections.Mapping)
    c = getattr(x, '_fields', None) is not None
    return a and b and c
