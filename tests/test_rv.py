import pytest
import raven
from raven.models.rv import RV, RVI, Ost, RVFile, isinstance_namedtuple
import datetime as dt
from collections import namedtuple
from .common import TESTDATA
from pathlib import Path


class TestRVFile:

    def test_simple_rv(self):
        fn = list(TESTDATA['raven-hmets'].glob('*.rvp'))[0]
        rvf = RVFile(fn)

        assert rvf.ext == 'rvp'
        assert rvf.stem == 'raven-hmets-salmon'
        assert not rvf.is_tpl

    def test_simple_tpl(self):
        fn = list(TESTDATA['ostrich-gr4j-cemaneige'].glob('*.rvp.tpl'))[0]
        rvf = RVFile(fn)

        assert rvf.ext == 'rvp'
        assert rvf.stem == 'raven-gr4j-salmon'
        assert rvf.is_tpl

    def test_ostIn(self):
        fn = list(TESTDATA['ostrich-gr4j-cemaneige'].glob('ostIn.txt'))[0]
        rvf = RVFile(fn)

        assert rvf.ext == 'txt'
        assert rvf.stem == 'ostIn'
        assert rvf.is_tpl

    def test_tags(self):
        rvp = list((Path(raven.__file__).parent / 'models' / 'raven-gr4j-cemaneige').glob("*.rvp"))[0]
        rvf = RVFile(rvp)

        assert isinstance(rvf.tags, list)
        assert 'params.GR4J_X3' in rvf.tags

    def test_fail(self):
        fn = Path(raven.__file__).parent
        with pytest.raises(ValueError):
            RVFile(fn)


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


class TestOst:
    def test_random(self):
        o = Ost()
        assert o.comment_random == '#'

        o = Ost(random_seed=0)
        assert o.comment_random == ''
        assert o.random_seed == 0


def test_isinstance_namedtuple():
    X = namedtuple('params', 'x1, x2, x3')
    x = X(1, 2, 3)
    assert isinstance_namedtuple(x)
    assert not isinstance_namedtuple([1, 2, 3])
