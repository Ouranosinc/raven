import pytest
from raven.models import RV
import datetime as dt
from collections import namedtuple

class TestRV:

    def test_end_date(self):
        rvi = RV(run_name='test',
                  start_date=dt.datetime(2000, 1, 1),
                  end_date=dt.datetime(2000, 1, 11),
                  )

        assert 10 == rvi.duration

        rvi.duration = 11
        assert dt.datetime(2000, 1, 12) == rvi.end_date

    def test_params(self):
        class MyClass(RV):
            params = namedtuple('p', 'x, y')

        rv = MyClass()
        rv['params'] = MyClass.params(1, 2)
        assert rv.params.x == 1

    def test_dict_interface(self):
        rv = RV(run_name='test')

        assert rv['run_name'] == rv.run_name

        with pytest.raises(AttributeError):
            rv['r'] = 6

    def test_flatten(self):
        class MyClass(RV):
            params = namedtuple('p', 'x, y')

        rv = MyClass(params=MyClass.params(1, 2))
        d = rv.flatten()
        assert 'x' in d

    def test_evaluation_metrics(self):
        rv = RV()
        rv.evaluation_metrics = 'LOG_NASH'

        with pytest.raises(ValueError):
            rv.evaluation_metrics = 'JIM'
