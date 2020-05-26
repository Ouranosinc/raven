import pytest
import raven
from raven.models.rv import (
    RV,
    RVI,
    RVT,
    RVC,
    Ost,
    RVFile,
    isinstance_namedtuple,
    RavenNcData,
    MonthlyAverage,
)
import datetime as dt
from collections import namedtuple
from .common import TESTDATA
from pathlib import Path
from io import StringIO


class TestRVFile:
    def test_simple_rv(self):
        fn = list(TESTDATA["raven-hmets"].glob("*.rvp"))[0]
        rvf = RVFile(fn)

        assert rvf.ext == "rvp"
        assert rvf.stem == "raven-hmets-salmon"
        assert not rvf.is_tpl

    def test_simple_tpl(self):
        fn = list(TESTDATA["ostrich-gr4j-cemaneige"].glob("*.rvp.tpl"))[0]
        rvf = RVFile(fn)

        assert rvf.ext == "rvp"
        assert rvf.stem == "raven-gr4j-salmon"
        assert rvf.is_tpl

    def test_ostIn(self):
        fn = list(TESTDATA["ostrich-gr4j-cemaneige"].glob("ostIn.txt"))[0]
        rvf = RVFile(fn)

        assert rvf.ext == "txt"
        assert rvf.stem == "ostIn"
        assert rvf.is_tpl

    def test_tags(self):
        rvp = list(
            (Path(raven.__file__).parent / "models" / "raven-gr4j-cemaneige").glob(
                "*.rvp"
            )
        )[0]
        rvf = RVFile(rvp)

        assert isinstance(rvf.tags, list)
        assert "params.GR4J_X3" in rvf.tags

    def test_fail(self):
        fn = Path(raven.__file__).parent
        with pytest.raises(ValueError):
            RVFile(fn)


class TestRV:
    def test_end_date(self):
        rvi = RVI(
            run_name="test",
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2000, 1, 11),
        )

        assert 10 == rvi.duration

        rvi.duration = 11
        assert dt.datetime(2000, 1, 12) == rvi.end_date

    def test_params(self):
        class RVP(RV):
            params = namedtuple("p", "x, y")

        rvp = RVP()
        rvp.params = RVP.params(1, 2)
        assert rvp.params.x == 1

    def test_dict_interface(self):
        rv = RV(run_name="test")

        assert rv["run_name"] == rv.run_name

        with pytest.raises(AttributeError):
            rv["r"] = 6

    def test_evaluation_metrics(self):
        rvi = RVI()
        rvi.evaluation_metrics = "LOG_NASH"

        with pytest.raises(ValueError):
            rvi.evaluation_metrics = "JIM"

    def test_update(self):
        rv = RV(a=None, b=None)
        rv.update({"a": 1, "b": 2})
        assert rv.a == 1

        rv.c = 1
        assert rv["c"] == 1

    def test_namedtuple(self):
        class Mod(RV):
            params = namedtuple("params", "x1, x2, x3")

        m = Mod(params=Mod.params(1, 2, 3))
        assert m.params.x1 == 1


def compare(a, b):
    """
    Compare two base strings, disregarding whitespace
    """
    import re

    return re.sub(r"\s*", "", a) == re.sub(r"\s*", "", b)


class TestRavenNcData:
    def test_simple(self):
        v = RavenNcData(
            var="tasmin",
            path="/path/tasmin.nc",
            var_name="tn",
            unit="deg_C",
            dimensions=["time", ],
        )
        tmp = str(v)

        assert compare(
            tmp,
            """:Data TEMP_MIN deg_C
                                  :ReadFromNetCDF
                                     :FileNameNC      /path/tasmin.nc
                                     :VarNameNC       tn
                                     :DimNamesNC      time
                                     :StationIdx      1
                                  :EndReadFromNetCDF
                               :EndData""",
        )

    def test_linear_transform(self):
        v = RavenNcData(
            var="tasmin",
            path="/path/tasmin.nc",
            var_name="tn",
            unit="deg_C",
            dimensions=["time", ],
            linear_transform=(24000.0, 0.0),
        )

        assert ":LinearTransform 24000.000000000000000 0.000000000000000" in str(v)


class TestMonthlyAve:
    def test_simple(self):
        ave = str(MonthlyAverage("Evaporation", range(12)))
        assert ave.startswith(":MonthlyAveEvaporation, 0, 1, 2")


class TestOst:
    def test_random(self):
        o = Ost()
        assert o.random_seed == ""

        o.random_seed = 0
        assert o.random_seed == "RandomSeed 0"


class TestRVI:
    def test_supress_output(self):
        rvi = RVI(suppress_output=True)
        assert rvi.suppress_output == ":SuppressOutput\n:DontWriteWatershedStorage"

        rvi = RVI(suppress_output=False)
        assert rvi.suppress_output == ""


rvc = TESTDATA["solution.rvc"].read_text()


class TestRVC:
    @classmethod
    def setup_class(self):
        self.rvc = StringIO(rvc)
        self.r = RVC()
        self.r.parse(self.rvc)

    def test_parse(self):
        assert self.r.hru_state.atmosphere == 821.98274
        assert self.r.basin_state.qout == [
            13.21660,
        ]
        assert self.r.basin_state.qoutlast == 13.29232

    def test_write(self):
        assert self.r.txt_hru_state.startswith("1,")
        assert self.r.txt_basin_state.strip().startswith(":BasinIndex 1,watershed")

    def test_format(self):
        rvc_template = Path(raven.models.__file__).parent / "global" / "global.rvc"
        params = dict(self.r.items())
        rvc_template.read_text().format(**params)


def test_isinstance_namedtuple():
    X = namedtuple("params", "x1, x2, x3")
    x = X(1, 2, 3)
    assert isinstance_namedtuple(x)
    assert not isinstance_namedtuple([1, 2, 3])
