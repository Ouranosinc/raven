import pytest
from . common import TESTDATA
from raven.models import Raven, GR4JCN, HMETS
import tempfile
import datetime as dt
import numpy as np


class TestRaven:

    def test_gr4j(self):
        rvs = TESTDATA['raven-gr4j-cemaneige-nc-rv']
        ts = TESTDATA['raven-gr4j-cemaneige-nc-ts']

        model = Raven(tempfile.mkdtemp())
        model.configure(rvs)
        model.run(ts)

    def test_mohyse(self):
        rvs = TESTDATA['raven-mohyse-rv']
        ts = list(TESTDATA['raven-mohyse-ts'])

        model = Raven(tempfile.mkdtemp())
        model.configure(rvs)
        model.run(ts)

    def test_hmets(self):
        rvs = TESTDATA['raven-hmets-rv']
        ts = TESTDATA['raven-hmets-ts']

        model = Raven(tempfile.mkdtemp())
        model.configure(rvs)
        model.run(ts)

    def test_hbvec(self):
        rvs = TESTDATA['raven-hbv-ec-rv']
        ts = TESTDATA['raven-hbv-ec-ts']

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
        model = GR4JCN()
        model.assign('run_name', 'test')
        assert model.rvi.run_name == 'test'

        model.assign('params', np.array([0.529, -3.396, 407.29, 1.072, 16.9, 0.947]))
        assert model.rvp.params.GR4J_X1 == 0.529

        model.assign('params', [0.529, -3.396, 407.29, 1.072, 16.9, 0.947])
        assert model.rvp.params.GR4J_X1 == 0.529

        model.assign('params', (0.529, -3.396, 407.29, 1.072, 16.9, 0.947))
        assert model.rvp.params.GR4J_X1 == 0.529

    def test_run(self):
        ts = TESTDATA['raven-gr4j-cemaneige-nc-ts']
        model = GR4JCemaneige()
        model(ts,
              start_date=dt.datetime(2000, 1, 1),
              end_date=dt.datetime(2002, 1, 1),
              area=4250.6,
              elevation=843.0,
              latitude=54.4848,
              longitude=-123.3659,
              params=(0.529, -3.396, 407.29, 1.072, 16.9, 0.947)
              )
        d = model.diagnostics
        np.testing.assert_almost_equal(d['DIAG_NASH_SUTCLIFFE'], -0.130614, 2)


class TestHMETS:

    def test_simple(self):
        ts = TESTDATA['raven-hmets-nc-ts']
        model = HMETS()
        params = (9.5019, 0.2774, 6.3942, 0.6884, 1.2875, 5.4134, 2.3641, 0.0973, 0.0464, 0.1998, 0.0222, -1.0919,
                  2.6851, 0.3740, 1.0000, 0.4739, 0.0114, 0.0243, 0.0069, 310.7211, 916.1947)

        model(ts,
              start_date=dt.datetime(2000, 1, 1),
              end_date=dt.datetime(2002, 1, 1),
              area=4250.6,
              elevation=843.0,
              latitude=54.4848,
              longitude=-123.3659,
              params=params
              )

        d = model.diagnostics
        np.testing.assert_almost_equal(d['DIAG_NASH_SUTCLIFFE'], -2.98165, 4)
