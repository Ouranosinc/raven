import pytest
from . common import TESTDATA
from raven.models import Raven, GR4JCemaneige, RVI, RV
import tempfile
import datetime as dt

class TestRaven:

    def test_gr4j(self):
        rvs = TESTDATA['raven-gr4j-cemaneige-nc-rv']
        ts = TESTDATA['raven-gr4j-cemaneige-nc-ts']

        model = Raven(tempfile.mkdtemp())
        model.configure(rvs)
        model.run([ts,], )

    def test_hmets(self):
        rvs = TESTDATA['raven-hmets-rv']
        ts = list(TESTDATA['raven-hmets-ts'])

        model = Raven(tempfile.mkdtemp())
        model.configure(rvs)
        model.run(ts)


class TestGR4JCemaneige:

    def test_simple(self):
        ts = TESTDATA['raven-gr4j-cemaneige-nc-ts']

        model = GR4JCemaneige(tempfile.mkdtemp())

        model.rvi.start_date = dt.datetime(2000, 1, 1)
        model.rvi.end_date   = dt.datetime(2002, 1, 1)
        model.rvi.run_name   = 'test'

        model.rvh.name      = 'Salmon'
        model.rvh.area      = '4250.6'
        model.rvh.elevation = '843.0'
        model.rvh.latitude  = 54.4848
        model.rvh.longitude = -123.3659

        rvp_set = [0.696, 0.7, 19.7, 2.09, 123.3, 0.25]
        # some params in Raven input files are derived from those 21 parameters
        # pdefaults.update({'GR4J_X1_hlf':            pdefaults['GR4J_X1']*1000./2.0})    --> x1 * 1000. / 2.0
        # pdefaults.update({'one_minus_CEMANEIGE_X2': 1.0 - pdefaults['CEMANEIGE_X2']})   --> 1.0 - x6
        x1       = rvp_set[0]
        x6       = rvp_set[5]
        rvp_set.append( x1 * 1000. / 2.0 )
        rvp_set.append( 1.0 - x6 )
        
        model.rvp.set(rvp_set)
        model.rvc.set(rvp_set)

        model.run([ts,])


    def test_tags(self):
        model = GR4JCemaneige(tempfile.mkdtemp())

        tags = model.tags
        assert 'run_name' in tags['rvi']


    def test_rvobjs(self):
        model = GR4JCemaneige(tempfile.mkdtemp())
        a = model.rvobjs
        a


    def test_assign(self):
        model = GR4JCemaneige(tempfile.mkdtemp())
        model.assign('run_name', 'test')
        assert model.rvi.run_name == 'test'
