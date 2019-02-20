from . common import TESTDATA
from raven.models import Raven, GR4JCN, HMETS, MOHYSE, HBVEC, GR4JCN_OST, HMETS_OST, MOHYSE_OST, HBVEC_OST
from raven.models import RavenMultiModel
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

        model([ts, ])

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

    def test_parallel(self):
        ts = TESTDATA['raven-gr4j-cemaneige-nc-ts']
        model = GR4JCN()
        model(ts,
              start_date=dt.datetime(2000, 1, 1),
              end_date=dt.datetime(2002, 1, 1),
              area=4250.6,
              elevation=843.0,
              latitude=54.4848,
              longitude=-123.3659,
              params=[(0.529, -3.396, 407.29, 1.072, 16.9, 0.947), (0.528, -3.4, 407.3, 1.07, 17, .95)]
              )

        assert len(model.diagnostics) == 2
        assert len(model.hydrograph) == 2


class TestGR4JCN_OST:
    def test_simple(self):
        ts = TESTDATA['ostrich-gr4j-cemaneige-nc-ts']
        model = GR4JCN_OST(test=True)
        params = (0.529, -3.396, 407.29, 1.072, 16.9, 0.053)
        low = (0.01, -15.0, 10.0, 0.0, 1.0, 0.0)
        high = (2.5, 10.0, 700.0, 7.0, 30.0, 1.0)

        model(ts,
              start_date=dt.datetime(1954, 1, 1),
              duration=208,
              area=4250.6,
              elevation=843.0,
              latitude=54.4848,
              longitude=-123.3659,
              params=params,
              lowerBounds=low,
              upperBounds=high,
              algorithm='DDS',
              random_seed='0',
              max_iterations=10,
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
        model = HMETS_OST(test=True)
        params = (9.5019, 0.2774, 6.3942, 0.6884, 1.2875, 5.4134, 2.3641, 0.0973, 0.0464, 0.1998, 0.0222, -1.0919,
                  2.6851, 0.3740, 1.0000, 0.4739, 0.0114, 0.0243, 0.0069, 310.7211, 916.1947)
        low = (0.3, 0.01, 0.5, 0.15, 0.0, 0.0, -2.0, 0.01, 0.0, 0.01, 0.005, -5.0, 0.0, 0.0, 0.0, 0.0,
               0.00001, 0.0, 0.00001, 0.0, 0.0)
        high = (20.0, 5.0, 13.0, 1.5, 20.0, 20.0, 3.0, 0.2, 0.1, 0.3, 0.1, 2.0, 5.0, 1.0, 3.0, 1.0,
                0.02, 0.1, 0.01, 0.5, 2.0)

        model(ts,
              start_date=dt.datetime(1954, 1, 1),
              duration=208,
              area=4250.6,
              elevation=843.0,
              latitude=54.4848,
              longitude=-123.3659,
              params=params,
              lowerBounds=low,
              upperBounds=high,
              algorithm='DDS',
              random_seed='0',
              max_iterations=10,
              )

        d = model.diagnostics
        np.testing.assert_almost_equal(d['DIAG_NASH_SUTCLIFFE'], -1.430270, 4)

        opt_para = model.calibrated_params
        opt_func = model.obj_func

        # # Random number seed: 123
        # # Budget:             50
        # # Algorithm:          DDS
        # # :StartDate          1954-01-01 00:00:00
        # # :Duration           20819
        # np.testing.assert_almost_equal( opt_para, [0.3243268,3.034247,407.2890,2.722774,12.18124,0.9468769], 4,
        #                                 err_msg='calibrated parameter set is not matching expected value')
        # np.testing.assert_almost_equal( opt_func, -0.5779910, 4,
        #                                 err_msg='calibrated NSE is not matching expected value')
        #
        # Random number seed: 123                         #
        # Budget:             10                          #      This is the setup used for testing:
        # Algorithm:          DDS                         #         shorter sim-period and lower budget
        # :StartDate          1954-01-01 00:00:00         #      First tested that example below matches
        # :Duration           208                         #
        np.testing.assert_almost_equal(opt_para, [1.777842E+01, 3.317211E+00, 5.727342E+00, 1.419491E+00,
                                                  1.382141E+01, 1.637954E+01, 7.166296E-01, 1.389346E-01,
                                                  2.620464E-02, 2.245525E-01, 2.839426E-02, -2.003810E+00,
                                                  9.479623E-01, 4.803857E-01, 2.524914E+00, 4.117232E-01,
                                                  1.950058E-02, 4.494123E-02, 1.405815E-03, 2.815803E-02,
                                                  1.007823E+00], 4,
                                       err_msg='calibrated parameter set is not matching expected value')
        np.testing.assert_almost_equal(opt_func, 1.430270, 4,
                                       err_msg='calibrated NSE is not matching expected value')

        # # Random number seed: 123                       #
        # # Budget:             50                        #      This is the setup in the Wiki:
        # # Algorithm:          DDS                       #      https://github.com/Ouranosinc/raven/wiki/
        # # :StartDate          1954-01-01 00:00:00       #      Technical-Notes#example-setups-for-hmets
        # # :Duration           20819                     #
        # np.testing.assert_almost_equal(opt_para, [5.008045E+00, 7.960246E-02, 4.332698E+00, 4.978125E-01,
        #                                           1.997029E+00, 6.269773E-01, 1.516961E+00, 8.180383E-02,
        #                                           6.730663E-02, 2.137822E-02, 2.097163E-02, 1.773348E+00,
        #                                           3.036039E-01, 1.928524E-02, 1.758471E+00, 8.942299E-01,
        #                                           8.741980E-03, 5.036474E-02, 9.465804E-03, 1.851839E-01,
        #                                           1.653934E-01, 2.624006E+00, 8.868485E-02, 9.259195E+01,
        #                                           8.269670E+01], 4,
        #                                err_msg='calibrated parameter set is not matching expected value')
        # np.testing.assert_almost_equal(opt_func, -6.350490E-01, 4,
        #                                err_msg='calibrated NSE is not matching expected value')


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
        np.testing.assert_almost_equal(d['DIAG_NASH_SUTCLIFFE'], 0.194612, 4)


class TestMOHYSE_OST():
    def test_simple(self):
        ts = TESTDATA['ostrich-mohyse-nc-ts']
        model = MOHYSE_OST(test=True)
        params = (1.0, 0.0468, 4.2952, 2.658, 0.4038, 0.0621, 0.0273, 0.0453)
        hrus = (0.9039, 5.6167)

        low_p = (0.01, 0.01, 0.01, -5.00, 0.01, 0.01, 0.01, 0.01)
        high_p = (20.0, 1.0, 20.0, 5.0, 0.5, 1.0, 1.0, 1.0)
        low_h = (0.01, 0.01)
        high_h = (15.0, 15.0)

        model(ts,
              start_date=dt.datetime(1954, 1, 1),
              duration=208,
              area=4250.6,
              elevation=843.0,
              latitude=54.4848,
              longitude=-123.3659,
              params=params,
              hrus=hrus,
              lowerBounds=low_p,
              upperBounds=high_p,
              hruslowerBounds=low_h,
              hrusupperBounds=high_h,
              algorithm='DDS',
              random_seed='0',
              max_iterations=10,
              )

        d = model.diagnostics
        np.testing.assert_almost_equal(d['DIAG_NASH_SUTCLIFFE'], 0.3826810, 4)

        opt_para = model.calibrated_params
        opt_func = model.obj_func

        # # Random number seed: 123
        # # Budget:             50
        # # Algorithm:          DDS
        # # :StartDate          1954-01-01 00:00:00
        # # :Duration           20819
        # np.testing.assert_almost_equal( opt_para, [0.3243268,3.034247,407.2890,2.722774,12.18124,0.9468769], 4,
        #                                 err_msg='calibrated parameter set is not matching expected value')
        # np.testing.assert_almost_equal( opt_func, -0.5779910, 4,
        #                                 err_msg='calibrated NSE is not matching expected value')
        #
        # Random number seed: 123                         #
        # Budget:             10                          #      This is the setup used for testing:
        # Algorithm:          DDS                         #         shorter sim-period and lower budget
        # :StartDate          1954-01-01 00:00:00         #      First tested that example below matches
        # :Duration           208                         #
        np.testing.assert_almost_equal(opt_para, [7.721801E+00, 8.551484E-01, 1.774571E+01, 1.627677E+00,
                                                  7.702450E-02, 9.409600E-01, 6.941596E-01, 8.207870E-01,
                                                  8.154455E+00, 1.018226E+01], 4,
                                       err_msg='calibrated parameter set is not matching expected value')
        np.testing.assert_almost_equal(opt_func, -0.3826810, 4,
                                       err_msg='calibrated NSE is not matching expected value')

        # # Random number seed: 123                       #
        # # Budget:             50                        #      This is the setup in the Wiki:
        # # Algorithm:          DDS                       #      https://github.com/Ouranosinc/raven/wiki/
        # # :StartDate          1954-01-01 00:00:00       #      Technical-Notes#example-setups-for-mohyse
        # # :Duration           20819                     #
        # np.testing.assert_almost_equal(opt_para, [1.517286E+01, 7.112556E-01, 1.981243E+01, -4.193046E+00,
        #                                           1.791486E-01, 9.774897E-01, 5.353541E-01, 6.686806E-01,
        #                                           1.040908E+01, 1.132304E+01, 8.831552E-02], 4,
        #                                err_msg='calibrated parameter set is not matching expected value')
        # np.testing.assert_almost_equal(opt_func, -0.3857010, 4,
        #                                err_msg='calibrated NSE is not matching expected value')


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
        np.testing.assert_almost_equal(d['DIAG_NASH_SUTCLIFFE'], -0.075407, 4)


class TestHBVEC_OST():
    def test_simple(self):
        ts = TESTDATA['ostrich-hbv-ec-nc-ts']
        model = HBVEC_OST(test=True)
        params = (0.05984519, 4.072232, 2.001574, 0.03473693, 0.09985144, 0.506052, 3.438486, 38.32455, 0.4606565,
                  0.06303738, 2.277781, 4.873686, 0.5718813, 0.04505643, 0.877607, 18.94145, 2.036937, 0.4452843,
                  0.6771759, 1.141608, 1.024278)

        low = (-3.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.01, 0.05, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01,
               0.0, 0.05, 0.8, 0.8)
        high = (3.0, 8.0, 8.0, 0.1, 1.0, 1.0, 7.0, 100.0, 1.0, 0.1, 6.0, 5.0, 5.0, 0.2, 1.0, 30.0, 3.0,
                2.0, 1.0, 1.5, 1.5)

        model(ts,
              start_date=dt.datetime(1954, 1, 1),
              duration=208,
              area=4250.6,
              elevation=843.0,
              latitude=54.4848,
              longitude=-123.3659,
              params=params,
              lowerBounds=low,
              upperBounds=high,
              algorithm='DDS',
              random_seed=0,
              max_iterations=10,
              )

        d = model.diagnostics
        np.testing.assert_almost_equal(d['DIAG_NASH_SUTCLIFFE'], -2.842420E-01, 4)

        opt_para = model.calibrated_params
        opt_func = model.obj_func

        # Random number seed: 123                         #
        # Budget:             10                          #      This is the setup used for testing:
        # Algorithm:          DDS                         #         shorter sim-period and lower budget
        # :StartDate          1954-01-01 00:00:00         #      First tested that example below matches
        # :Duration           208                         #
        np.testing.assert_almost_equal(opt_para, [-8.317931E-01, 4.072232E+00, 2.001574E+00, 5.736299E-03,
                                                  9.985144E-02, 4.422529E-01, 3.438486E+00, 8.055843E+01,
                                                  4.440133E-01, 8.451082E-02, 2.814201E+00, 7.327970E-01,
                                                  1.119773E+00, 1.161223E-03, 4.597179E-01, 1.545857E+01,
                                                  1.223865E+00, 4.452843E-01, 9.492006E-01, 9.948123E-01,
                                                  1.110682E+00], 4,
                                       err_msg='calibrated parameter set is not matching expected value')
        np.testing.assert_almost_equal(opt_func, 2.842420E-01, 4,
                                       err_msg='calibrated NSE is not matching expected value')

        # # Random number seed: 123                       #
        # # Budget:             50                        #      This is the setup in the Wiki:
        # # Algorithm:          DDS                       #      https://github.com/Ouranosinc/raven/wiki/
        # # :StartDate          1954-01-01 00:00:00       #      Technical-Notes#example-setups-for-environment-
        # # :Duration           20819                     #
        # np.testing.assert_almost_equal(opt_para, [5.984519E-02, 4.072232E+00, 2.001574E+00, 3.473693E-02,
        #                                           9.985144E-02, 5.060520E-01, 2.944343E+00, 3.832455E+01,
        #                                           4.606565E-01, 6.303738E-02, 2.277781E+00, 4.873686E+00,
        #                                           5.718813E-01, 4.505643E-02, 8.776511E-01, 1.894145E+01,
        #                                           2.036937E+00, 4.452843E-01, 6.771759E-01, 1.206053E+00,
        #                                           1.024278E+00], 4,
        #                                err_msg='calibrated parameter set is not matching expected value')
        # np.testing.assert_almost_equal(opt_func, -6.034670E-01, 4,
        #                                err_msg='calibrated NSE is not matching expected value')


class TestRavenMultiModel:

    def test_simple(self):
        ts = TESTDATA['raven-hmets-nc-ts']
        model = RavenMultiModel(models=['gr4jcn', 'hmets'])
        gr4jcn = (0.529, -3.396, 407.29, 1.072, 16.9, 0.947)
        hmets = (9.5019, 0.2774, 6.3942, 0.6884, 1.2875, 5.4134, 2.3641, 0.0973, 0.0464, 0.1998, 0.0222, -1.0919,
                 2.6851, 0.3740, 1.0000, 0.4739, 0.0114, 0.0243, 0.0069, 310.7211, 916.1947)

        model(ts,
              start_date=dt.datetime(2000, 1, 1),
              end_date=dt.datetime(2002, 1, 1),
              area=4250.6,
              elevation=843.0,
              latitude=54.4848,
              longitude=-123.3659,
              gr4jcn=gr4jcn,
              hmets=hmets,
              )

        assert len(model.q_sim) == 2
