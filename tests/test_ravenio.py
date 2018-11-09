import pytest
from raven.processes import ravenio
from pathlib import Path
from .common import TESTDATA
import datetime as dt
import os
import tempfile
from collections import OrderedDict as odict
import glob

"""
class TestStartDate:

    def test_grj4_cemaneige(self):
        start, end = ravenio.start_end_date(TESTDATA['gr4j-cemaneige'].values())
        assert start == dt.datetime(2000, 7, 1)
        assert end == dt.datetime(2002, 7, 1)

    def test_raven_gr4j_cemaneige(self):
        start, end = ravenio.start_end_date([TESTDATA['raven-gr4j-cemaneige-nc-ts'],])
        assert start == dt.datetime(1954, 1, 1)
        assert end == dt.datetime(2010, 12, 31)

pytest.mark.skip
class TestSetupModel:
    gr4j = odict(
        rvi=dict(run_name='run1', Start_Date='2000-01-01 00:00:00', End_Date=None, Duration=10, TimeStep=1.0, EvaluationMetrics='NASH_SUTCLIFFE RMSE'),
        rvp=odict(GR4J_X1=.7, GR4J_X2=.7, GR4J_X3=19., GR4J_X4=2.09, AvgAnnualSnow=123, AirSnowCoeff=.75),
        rvc=odict(SOIL_0=1, SOIL_1=2),
        rvh=dict(NAME='Test', AREA=45, ELEVATION=3, LATITUDE=45, LONGITUDE=-154),
        rvt=dict(pr='data/data.nc', prsn='data/data.nc', tasmin='data/data.nc', tasmax='data/data.nc',
                 evspsbl='data/data.nc', water_volume_transport_in_river_channel='data/data.nc',
                 pr_var='rain', prsn_var='snow', tasmin_var='tmin', tasmax_var='tmax', evspsbl_var='pet',
                 water_volume_transport_in_river_channel_var='qobs')
    )

    def test_gr4j(self):
        name = 'raven-gr4j-cemaneige'
        outpath = tempfile.mkdtemp()
        ravenio.setup_model(name, outpath, self.gr4j)

        # Make sure there are no {} in the filled templates
        for fn in glob.glob(os.path.join(outpath, 'model', '*.rv?')):
            with open(fn) as f:
                txt = f.read()
                assert '{' not in txt
                assert '}' not in txt

@pytest.mark.skip
class TestReadDiagnostics:
    def test_simple(self):
        import StringIO
        diag = StringIO.StringIO()
        diag.write('observed data series,filename,DIAG_NASH_SUTCLIFFE,DIAG_RMSE,\n'
                   'HYDROGRAPH,data_obs/Salmon-River-Near-Prince-George_meteo_daily.nc,-0.746096,62.1502,')
        diag.seek(0)

        out = ravenio.read_diagnostics(diag)
        assert out['DIAG_RMSE'] == 62.1502
"""

class TestParseConfiguration:
    def test_simple(self):
        p = Path(TESTDATA['raven-hmets'])
        rvi = list(p.glob("*.rvi"))[0]

        out = ravenio.parse_configuration(rvi)
        assert out['Duration'] == '2081'

