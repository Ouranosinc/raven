from . common import TESTDATA
from raven.models import Raven, Ostrich
import tempfile


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
        assert len(model.calibrated_params) == 6
        assert isinstance(model.obj_func, float)
