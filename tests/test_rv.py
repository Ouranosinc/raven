import pytest
from raven.models.rv import RV, RVI, isinstance_namedtuple
import datetime as dt
from collections import namedtuple


class TestRV:

    def test_end_date(self):
        rvi = RVI(run_name='test',
                  start_date=dt.datetime(2000, 1, 1),
                  end_date=dt.datetime(2000, 1, 11),
                  )

        assert 10 == rvi.duration

        rvi.duration = 11
        assert dt.datetime(2000, 1, 12) == rvi.end_date

    def test_params(self):
        class RVP(RV):
            params = namedtuple('p', 'x, y')

        rvp = RVP()
        rvp.params = RVP.params(1, 2)
        assert rvp.params.x == 1

        d = rvp.to_dict()
        assert 'x' in d

    def test_dict_interface(self):
        rv = RV(run_name='test')

        assert rv['run_name'] == rv.run_name

        with pytest.raises(AttributeError):
            rv['r'] = 6

    def test_evaluation_metrics(self):
        rvi = RVI()
        rvi.evaluation_metrics = 'LOG_NASH'

        with pytest.raises(ValueError):
            rvi.evaluation_metrics = 'JIM'

    def test_update(self):
        rv = RV(a=None, b=None)
        rv.update({'a': 1, 'b': 2})
        assert rv.a == 1

        rv.c = 1
        assert rv['c'] == 1

    def test_namedtuple(self):
        class Mod(RV):
            params = namedtuple('params', 'x1, x2, x3')

        m = Mod(params=Mod.params(1, 2, 3))
        assert m.params.x1 == 1


def test_isinstance_namedtuple():
    X = namedtuple('params', 'x1, x2, x3')
    x = X(1, 2, 3)
    assert isinstance_namedtuple(x)
    assert not isinstance_namedtuple([1, 2, 3])
