import six
import datetime as dt
from collections import OrderedDict





class RV(object):
    _keys = None
    _magic_key_name = ''

    def __init__(self, **kwargs):
        if self._keys is None:
            self._keys = self.__slots__ = tuple(kwargs.keys())

        for attribute in self._keys:
            setattr(self, attribute, kwargs.pop(attribute))

    def keys(self):
        """Return a list of keys."""
        return self._keys

    def values(self):
        """Return a list of values."""
        return [getattr(self, key) for key in self._keys]

    def set(self, x):
        """Set the values from a list, in the same order as this object's keys."""
        for key, val in zip(self._keys, x):
            setattr(self, key, val)

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        if key == self._magic_key_name:
            self.set(value)
        elif key in self._keys:
            setattr(self, key, value)
        else:
            raise AttributeError("Cannot create attribute {}".format(key))

    def __len__(self):
        return len(self._keys)

    def __iter__(self):
        return iter(self._keys)

    def items(self):
        """Return a generator of (key, value) pairs."""
        for attribute in self._keys:
            yield attribute, getattr(self, attribute)


class RVI(RV):
    _keys = ('run_name', 'start_date', 'end_date', 'duration', 'time_step', 'evaluation_metrics',)

    def __init__(self,
                 run_name='raven-sim',
                 start_date=None,
                 end_date=None,
                 duration=1,
                 time_step=1.0,
                 evaluation_metrics='NASH_SUTCLIFFE RMSE',
                 **kwargs):

        self._run_name = run_name
        self._start_date = start_date
        self._end_date = end_date
        self._duration = duration
        self._time_step = time_step
        self._evaluation_metrics = evaluation_metrics

        for key, val in kwargs.items():
            setattr(self, key, val)

        self._update_duration()
        self._update_end_date()

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


class RVP(RV):
    _magic_key_name = 'params'


class RVC(RV):
    _magic_key_name = 'init'


class RVH(RV):
    _magic_key_name = 'hru'
