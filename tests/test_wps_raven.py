import pytest
import os
import glob
#import pathlib
from pywps import Service
from pywps.tests import assert_response_success

from . common import client_for, TESTDATA, CFG_FILE
from raven.processes import RavenProcess

cf = ['rvi', 'rvp', 'rvc', 'rvh', 'rvt']
pattern = TESTDATA['raven-gr4j-cemaneige-nc-rv']
rvfiles = {os.path.splitext(fn)[1][1:]: os.path.join(os.path.split(pattern)[0], fn) for fn in glob.glob(pattern)}
for f in rvfiles.values():
    if not os.path.exists(f):
        raise ValueError("{} does not exist".format(f))

class TestRavenProcess:

    def test_gr4j_saumon_nc(self):
        client = client_for(Service(processes=[RavenProcess(), ], cfgfiles=CFG_FILE))

        datainputs = ("nc=files@xlink:href=file://{fn};" + \
                     ';'.join(["conf=files@xlink:href=file://{%s}"%key for key
                                                                     in cf])) \
            .format(fn=TESTDATA['raven-gr4j-cemaneige-nc-ts'], **rvfiles)

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='raven',
            datainputs=datainputs)
        print(resp.response[0])
        assert_response_success(resp)



"rvi=files@xlink:href=file://{rvi};" \
                     "rvp=files@xlink:href=file://{rvp};" \
                     "rvc=files@xlink:href=file://{rvc};" \
                     "rvh=files@xlink:href=file://{rvh};" \
                     "rvt=files@xlink:href=file://{rvt};" \
