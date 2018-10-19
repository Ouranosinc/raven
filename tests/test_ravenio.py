import pytest
from raven.processes import ravenio
from raven.models import raven_templates

from .common import TESTDATA
import datetime as dt
import os
import tempfile
from collections import OrderedDict as odict
import glob


class TestStartDate:

    def test_grj4_cemaneige(self):
        start, end = ravenio.start_end_date(TESTDATA['gr4j-cemaneige'].values())
        assert start == dt.datetime(2000, 7, 1)
        assert end == dt.datetime(2002, 7, 1)

    def test_raven_gr4j_cemaneige(self):
        start, end = ravenio.start_end_date([TESTDATA['raven-gr4j-cemaneige'],])
        assert start == dt.datetime(1954, 1, 1)
        assert end == dt.datetime(2010, 12, 31)

class TestSetupModel:
    gr4j = odict(
        rvi=dict(Start_Date='2000-01-01 00:00:00', End_Date=None, Duration=10, TimeStep=1.0, EvaluationMetrics='NASH_SUTCLIFFE RMSE'),
        rvp=odict(GR4J_X1=.7, GR4J_X2=.7, GR4J_X3=19., GR4J_X4=2.09, AvgAnnualSnow=123, AirSnowCoeff=.75),
        rvc=odict(SOIL_0=1, SOIL_1=2),
        rvh=dict(NAME='Test', AREA=45, ELEVATION=3, LATITUDE=45, LONGITUDE=-154),
        rvt=dict(RAIN='data/data.nc', SNOW='data/data.nc', TASMIN='data/data.nc', TASMAX='data/data.nc', PET='data/data.nc', QOBS='data/data.nc')
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


