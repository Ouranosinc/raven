import pytest
import raven
from raven.models.rv import RV, RVI, RVT, RVC, Ost, RVFile, isinstance_namedtuple, RavenNcData, MonthlyAverage
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


def compare(a, b):
    """
    Compare two base strings, disregarding whitespace
    """
    import re
    return re.sub(r"\s*", "", a) == re.sub(r"\s*", "", b)


class TestRavenNcData:

    def test_simple(self):
        v = RavenNcData(var='tasmin', path='/path/tasmin.nc', var_name='tn', unit='deg_C', dimensions=['time', ])
        tmp = str(v)

        assert compare(tmp, """:Data TEMP_MIN deg_C
                                  :ReadFromNetCDF
                                     :FileNameNC      /path/tasmin.nc
                                     :VarNameNC       tn
                                     :DimNamesNC      time
                                     :StationIdx      1
                                  :EndReadFromNetCDF
                               :EndData""")

    def test_linear_transform(self):
        v = RavenNcData(var='tasmin', path='/path/tasmin.nc', var_name='tn', unit='deg_C', dimensions=['time', ],
                        linear_transform=(24000.0, 0.0))

        assert ":LinearTransform 24000.000000000000000 0.000000000000000" in str(v)


class TestMonthlyAve:
    def test_simple(self):
        ave = str(MonthlyAverage('Evaporation', range(12)))
        assert ave.startswith(':MonthlyAveEvaporation, 0, 1, 2')


class TestOst:
    def test_random(self):
        o = Ost()
        assert o.random_seed == ''

        o.random_seed = 0
        assert o.random_seed == 'RandomSeed 0'


class TestRVI:
    def test_supress_output(self):
        rvi = RVI(suppress_output=True)
        assert rvi.suppress_output == ":SuppressOutput\n:DontWriteWatershedStorage"

        rvi = RVI(suppress_output=False)
        assert rvi.suppress_output == ""


rvc = \
""":TimeStamp 2002-01-01 00:00:00.00
:HRUStateVariableTable
 :Attributes,SURFACE_WATER,ATMOSPHERE,ATMOS_PRECIP,PONDED_WATER,SOIL[0],SOIL[1],SOIL[2],SOIL[3],SNOW_TEMP,SNOW,SNOW_COVER,AET,CONVOLUTION[0],CONVOLUTION[1],CONV_STOR[0],CONV_STOR[1],CONV_STOR[2],CONV_STOR[3],CONV_STOR[4],CONV_STOR[5],CONV_STOR[6],CONV_STOR[7],CONV_STOR[8],CONV_STOR[9],CONV_STOR[10],CONV_STOR[11],CONV_STOR[12],CONV_STOR[13],CONV_STOR[14],CONV_STOR[15],CONV_STOR[16],CONV_STOR[17],CONV_STOR[18],CONV_STOR[19],CONV_STOR[20],CONV_STOR[21],CONV_STOR[22],CONV_STOR[23],CONV_STOR[24],CONV_STOR[25],CONV_STOR[26],CONV_STOR[27],CONV_STOR[28],CONV_STOR[29],CONV_STOR[30],CONV_STOR[31],CONV_STOR[32],CONV_STOR[33],CONV_STOR[34],CONV_STOR[35],CONV_STOR[36],CONV_STOR[37],CONV_STOR[38],CONV_STOR[39],CONV_STOR[40],CONV_STOR[41],CONV_STOR[42],CONV_STOR[43],CONV_STOR[44],CONV_STOR[45],CONV_STOR[46],CONV_STOR[47],CONV_STOR[48],CONV_STOR[49],CONV_STOR[50],CONV_STOR[51],CONV_STOR[52],CONV_STOR[53],CONV_STOR[54],CONV_STOR[55],CONV_STOR[56],CONV_STOR[57],CONV_STOR[58],CONV_STOR[59],CONV_STOR[60],CONV_STOR[61],CONV_STOR[62],CONV_STOR[63],CONV_STOR[64],CONV_STOR[65],CONV_STOR[66],CONV_STOR[67],CONV_STOR[68],CONV_STOR[69],CONV_STOR[70],CONV_STOR[71],CONV_STOR[72],CONV_STOR[73],CONV_STOR[74],CONV_STOR[75],CONV_STOR[76],CONV_STOR[77],CONV_STOR[78],CONV_STOR[79],CONV_STOR[80],CONV_STOR[81],CONV_STOR[82],CONV_STOR[83],CONV_STOR[84],CONV_STOR[85],CONV_STOR[86],CONV_STOR[87],CONV_STOR[88],CONV_STOR[89],CONV_STOR[90],CONV_STOR[91],CONV_STOR[92],CONV_STOR[93],CONV_STOR[94],CONV_STOR[95],CONV_STOR[96],CONV_STOR[97],CONV_STOR[98],CONV_STOR[99]
 :Units,mm,mm,mm,mm,mm,mm,mm,mm,C,mm,0-1,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm
 1,0.00000,821.98274,-1233.16000,0.00000,276.27919,124.08013,0.00000,47.52822,-8.71223,66.93150,1.00000,0.00000,0.02884,0.01171,-299.21993,0.02884,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,-33.24666,0.01165,0.00007,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000
:EndHRUStateVariableTable
:BasinStateVariables
 :BasinIndex 1,watershed
   :ChannelStorage, 0.00000
   :RivuletStorage, 570957.08719
   :Qout,1,13.21660,13.29232
   :Qlat,3,13.21660,13.29232,13.36898,13.21660
   :Qin ,20,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000
:EndBasinStateVariables
"""

class TestRVC:
    def test_parse(self):
        from io import StringIO
        solution = StringIO(rvc)
        r = RVC()
        r.parse(solution)

        assert r._hru_state[1].atmosphere == 821.98274
        assert r._basin_state[1].qout == [13.21660, 13.29232]


def test_isinstance_namedtuple():
    X = namedtuple('params', 'x1, x2, x3')
    x = X(1, 2, 3)
    assert isinstance_namedtuple(x)
    assert not isinstance_namedtuple([1, 2, 3])
