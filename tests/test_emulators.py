from . common import TESTDATA, make_bnds
from raven.models import Raven, GR4JCN, HMETS, MOHYSE, HBVEC, GR4JCN_OST, HMETS_OST, MOHYSE_OST, HBVEC_OST
import tempfile
import datetime as dt
import numpy as np


class TestGR4JCN:

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

        model.rvp.params = model.params(0.529, -3.396, 407.29, 1.072, 16.9, 0.947)

        model.run([ts, ])

        d = model.diagnostics
        # yields NSE=0.5112 for full period 1954-2010
        np.testing.assert_almost_equal(d['DIAG_NASH_SUTCLIFFE'], -0.130614, 2)

        hds = model.q_sim
        assert hds.attrs['long_name'] == 'Simulated outflows'

    def test_tags(self):
        model = GR4JCN(tempfile.mkdtemp())

        tags = model.tags
        assert 'run_name' in tags

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
        model = GR4JCN()
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

    def test_overwrite(self):
        ts = TESTDATA['raven-gr4j-cemaneige-nc-ts']
        model = GR4JCN()
        model(ts,
              start_date=dt.datetime(2000, 1, 1),
              end_date=dt.datetime(2002, 1, 1),
              area=4250.6,
              elevation=843.0,
              latitude=54.4848,
              longitude=-123.3659,
              params=(0.529, -3.396, 407.29, 1.072, 16.9, 0.947)
              )

        qsim1 = model.q_sim.copy(deep=True)
        m1 = qsim1.mean()

        model(ts, params=(0.528, -3.4, 407.3, 1.07, 17, .95), overwrite=True)

        qsim2 = model.q_sim.copy(deep=True)
        m2 = qsim2.mean()
        assert m1 != m2

        np.testing.assert_almost_equal(m1, m2, 1)

        d = model.diagnostics
        np.testing.assert_almost_equal(d['DIAG_NASH_SUTCLIFFE'], -0.130614, 2)

    def test_version(self):
        model = Raven()
        assert model.version == '2.9'

        model = GR4JCN()
        assert model.version == '2.9'


class TestGR4JCN_OST:
    def test_simple(self):
        ts = TESTDATA['ostrich-gr4j-cemaneige-nc-ts']
        model = GR4JCN_OST()
        low = (0.1, -5.0, 100.0, 1.0, 10.0, 0.1)
        high = (0.9, 0.0, 500.0, 1.1, 20.0, 1.0)

        model(ts,
              start_date=dt.datetime(1954, 1, 1),
              duration=208,
              area=4250.6,
              elevation=843.0,
              latitude=54.4848,
              longitude=-123.3659,
              lowerBounds=low,
              upperBounds=high,
              algorithm='DDS',
              MaxEvals=10,
              )

        d = model.diagnostics
        np.testing.assert_almost_equal(d['DIAG_NASH_SUTCLIFFE'], 0.5078130, 4)

        # Random number seed: 123
        # Budget:             10
        # Algorithm:          DDS
        # :StartDate          1954-01-01 00:00:00
        # :Duration           208
        opt_para = model.calibrated_params
        opt_func = model.obj_func
        np.testing.assert_almost_equal(opt_para, [2.424726, 3.758972, 204.3856, 5.866946, 16.60408, 0.3728098], 4,
                                       err_msg='calibrated parameter set is not matching expected value')
        np.testing.assert_almost_equal(opt_func, -0.5078130, 4,
                                       err_msg='calibrated NSE is not matching expected value')

        # # Random number seed: 123
        # # Budget:             50
        # # Algorithm:          DDS
        # # :StartDate          1954-01-01 00:00:00
        # # :Duration           20819
        # np.testing.assert_almost_equal( opt_para, [0.3243268,3.034247,407.2890,2.722774,12.18124,0.9468769], 4,
        #                                 err_msg='calibrated parameter set is not matching expected value')
        # np.testing.assert_almost_equal( opt_func, -0.5779910, 4,
        #                                 err_msg='calibrated NSE is not matching expected value')


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


class TestHMETS_OST:

    def test_simple(self):
        ts = TESTDATA['raven-hmets-nc-ts']
        model = HMETS_OST()
        params = (9.5019, 0.2774, 6.3942, 0.6884, 1.2875, 5.4134, 2.3641, 0.0973, 0.0464, 0.1998, 0.0222, -1.0919,
                  2.6851, 0.3740, 1.0000, 0.4739, 0.0114, 0.0243, 0.0069, 310.7211, 916.1947)
        low, high = make_bnds(params, .25)

        model(ts,
              start_date=dt.datetime(1954, 1, 1),
              duration=208,
              area=4250.6,
              elevation=843.0,
              latitude=54.4848,
              longitude=-123.3659,
              lowerBounds=low,
              upperBounds=high,
              algorithm='DDS',
              MaxEvals=10
              )


class TestMOHYSE:

    def test_simple(self):
        ts = TESTDATA['raven-mohyse-nc-ts']
        model = MOHYSE()
        params = (1.0, 0.0468, 4.2952, 2.658, 0.4038, 0.0621, 0.0273, 0.0453)
        hrus = (0.9039, 5.6167)

        model(ts,
              start_date=dt.datetime(2000, 1, 1),
              end_date=dt.datetime(2002, 1, 1),
              area=4250.6,
              elevation=843.0,
              latitude=54.4848,
              longitude=-123.3659,
              params=params,
              hrus=hrus,
              )

        d = model.diagnostics
        # TODO: Correct expected NSE
        np.testing.assert_almost_equal(d['DIAG_NASH_SUTCLIFFE'], 0.194612, 4)


class TestMOHYSE_OST():
    def test_simple(self):
        ts = TESTDATA['ostrich-mohyse-nc-ts']
        model = MOHYSE_OST()
        params = (1.0, 0.0468, 4.2952, 2.658, 0.4038, 0.0621, 0.0273, 0.0453)
        hrus = (0.9039, 5.6167)

        low_p, high_p = make_bnds(params, .25)
        low_h, high_h = make_bnds(hrus, .5)

        model(ts,
              start_date=dt.datetime(1954, 1, 1),
              duration=208,
              area=4250.6,
              elevation=843.0,
              latitude=54.4848,
              longitude=-123.3659,
              lowerBounds=low_p,
              upperBounds=high_p,
              hruslowerBounds=low_h,
              hrusupperBounds=high_h,
              MaxEvals=10
              )

        d = model.diagnostics
        np.testing.assert_almost_equal(d['DIAG_NASH_SUTCLIFFE'], 0.5078130, 4)

        opt_para = model.calibrated_params
        opt_func = model.obj_func
        np.testing.assert_almost_equal(opt_para, [1.0, 0.0468, 4.2952, 2.658, 0.4038, 0.0621, 0.0273, 0.0453], 4,
                                       err_msg='calibrated parameter set is not matching expected value')
        np.testing.assert_almost_equal(opt_func, -0.5078130, 4,
                                       err_msg='calibrated NSE is not matching expected value')


class TestHBVEC:

    def test_simple(self):
        ts = TESTDATA['raven-hbv-ec-nc-ts']
        model = HBVEC()
        params = (0.05984519, 4.072232, 2.001574, 0.03473693, 0.09985144, 0.506052, 3.438486, 38.32455, 0.4606565,
                  0.06303738, 2.277781, 4.873686, 0.5718813, 0.04505643, 0.877607, 18.94145, 2.036937, 0.4452843,
                  0.6771759, 1.141608, 1.024278)

        model(ts,
              start_date=dt.datetime(2000, 1, 1),
              end_date=dt.datetime(2002, 1, 1),
              area=4250.6,
              elevation=843.0,
              latitude=54.4848,
              longitude=-123.3659,
              params=params,
              )

        d = model.diagnostics
        # TODO: Correct expected NSE
        np.testing.assert_almost_equal(d['DIAG_NASH_SUTCLIFFE'], -0.075407, 4)
