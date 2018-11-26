import pytest
from raven.models import RV, RVI, RVP
import datetime as dt


class TestRV:

    def test_rvi(self):
        rvi = RVI(run_name='test',
                  start_date=dt.datetime(2000, 1, 1),
                  end_date=dt.datetime(2000, 1, 11),
                  )

        assert 10 == rvi.duration

        rvi.duration = 11
        assert dt.datetime(2000, 1, 12) == rvi.end_date

    def test_rvp(self):
        rvp = RV(GR4J_X1=1, GR4J_X2=2, GR4J_X3=3, GR4J_X4=4, CEMANEIGE_X1=5, CEMANEIGE_X2=6)
        assert list(range(1, 7)) == rvp.values()

    def test_dict_interface(self):
        rvp = RV(a=1, b=2)
        d = dict(c=3)
        d.update(rvp.items())
        assert 'a' in d

        with pytest.raises(AttributeError):
            rvp['r'] = 6

    def test_magic(self):
        rvp = RVP(a=1, b=2)
        rvp['params'] = [3, 4]

        assert 3 == rvp.a
