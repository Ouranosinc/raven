from . common import TESTDATA
from raven.models import Raven, Ostrich
import tempfile
import numpy as np


class TestRaven:

    def test_gr4j(self):
        rvs = TESTDATA['raven-gr4j-cemaneige-nc-rv']
        ts = TESTDATA['raven-gr4j-cemaneige-nc-ts']

        model = Raven()
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


class TestOstrich:

    def test_gr4j_with_no_tags(self):
        ts = TESTDATA['raven-gr4j-cemaneige-nc-ts']
        ost = TESTDATA['ostrich-gr4j-cemaneige-rv']

        model = Ostrich()
        model.configure(ost)
        print(model.exec_path)
        model.run(ts)

        opt_para = model.calibrated_params
        opt_func = model.obj_func

        assert len(opt_para) == 6
        assert isinstance(opt_func, float)

        # Random number seed: 123
        # Budget:             10
        # Algorithm:          DDS
        # :StartDate          1954-01-01 00:00:00
        # :Duration           208
        np.testing.assert_almost_equal(opt_para, [2.424726, 3.758972, 204.3856, 5.866946, 16.60408, 0.3728098], 4,
                                       err_msg='calibrated parameter set is not matching expected value')
        np.testing.assert_almost_equal(opt_func, -0.5078130, 4,
                                       err_msg='calibrated NSE is not matching expected value')

        # # Random number seed: 123
        # # Budget:             50
        # # Algorithm:          DDS
        # # :StartDate          1954-01-01 00:00:00
        # # :Duration           20819
        # np.testing.assert_almost_equal(opt_para, [0.3243268,3.034247,407.2890,2.722774,12.18124,0.9468769], 4,
        #                                err_msg='calibrated parameter set is not matching expected value')
        # np.testing.assert_almost_equal(opt_func, -0.5779910, 4,
        #                                err_msg='calibrated NSE is not matching expected value')
