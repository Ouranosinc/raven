import pytest
from . common import TESTDATA
from raven.models import Raven
import tempfile

class TestRaven:

    def test_init(self):
        rvs = TESTDATA['raven-gr4j-cemaneige-nc-rv']
        ts = TESTDATA['raven-gr4j-cemaneige-nc-ts']

        r = Raven(tempfile.mkdtemp())
        r.configure(**{f.suffix[1:]: f for f in rvs})
        r.run([ts,], )
