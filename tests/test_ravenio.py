import pytest
from raven.processes import ravenio
from .common import TESTDATA
import datetime as dt

class TestStartDate:

    def test_grj4_cemaneige(self):
        start = ravenio.start_date(TESTDATA['gr4j-cemaneige'].values())
        assert start == dt.datetime(2000, 7, 1)


    def test_raven_gr4j_cemaneige(self):
        start = ravenio.start_date([TESTDATA['raven-gr4j-cemaneige'],])
        assert start == dt.datetime(1954, 1, 1)
