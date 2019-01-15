import pytest
from . common import TESTDATA
from raven.models import Raven, GR4JCN, RVI, RV
import tempfile
import datetime as dt
import numpy as np


class TestRaven:

    def test_gr4j(self):
        rvs = TESTDATA['raven-gr4j-cemaneige-nc-rv']
        ts = TESTDATA['raven-gr4j-cemaneige-nc-ts']

        model = Raven(tempfile.mkdtemp())
        model.configure(rvs)
        model.run([ts, ], )

    def test_mohyse(self):
        rvs = TESTDATA['raven-mohyse-rv']
        ts = list(TESTDATA['raven-mohyse-ts'])

        model = Raven(tempfile.mkdtemp())
        model.configure(rvs)
        model.run(ts)

    def test_hmets(self):
        rvs = TESTDATA['raven-hmets-rv']
        ts = list(TESTDATA['raven-hmets-ts'])

        model = Raven(tempfile.mkdtemp())
        model.configure(rvs)
        model.run(ts)

    def test_hbvec(self):
        rvs = TESTDATA['raven-hbv-ec-rv']
        ts = list(TESTDATA['raven-hbv-ec-ts'])

        model = Raven(tempfile.mkdtemp())
        model.configure(rvs)
        model.run(ts)


class TestGR4JCemaneige:

    def test_simple(self):
        ts = TESTDATA['raven-gr4j-cemaneige-nc-ts']

        model = GR4JCN(tempfile.mkdtemp())

        model.rvi.start_date = dt.datetime(2000, 1, 1)
        model.rvi.end_date = dt.datetime(2002, 1, 1)
        model.rvi.run_name = 'test'

        model.rvh.name = 'Salmon'
        model.rvh.area = '4250.6'
        model.rvh.elevation = '843.0'
        model.rvh.latitude = 54.4848
        model.rvh.longitude = -123.3659

        model.rvp.params = model.RVP.params(0.529, -3.396, 407.29, 1.072, 16.9, 0.947)

        model.run([ts, ])

        d = model.diagnostics
        # yields NSE=0.5112 for full period 1954-2010
        np.testing.assert_almost_equal(d['DIAG_NASH_SUTCLIFFE'], -0.130614, 2)

        hds = model.hydrograph
        assert 'q_sim' in hds.data_vars

    def test_tags(self):
        model = GR4JCN(tempfile.mkdtemp())

        tags = model.tags
        assert 'run_name' in tags['rvi']

    def test_rvobjs(self):
        model = GR4JCN(tempfile.mkdtemp())
        a = model.rvobjs
        assert a

    def test_assign(self):
        model = GR4JCN(tempfile.mkdtemp())
        model.assign('run_name', 'test')
        assert model.rvi.run_name == 'test'
