import datetime as dt
import os
import tempfile
import numpy as np
import pytest

from raven.utilities.data_assimilation import assimilateQobsSingleDay
from raven.models import GR4JCN
from .common import TESTDATA

class TestAssimilationGR4JCN:
    def test_simple(self):
        ts = TESTDATA["raven-gr4j-cemaneige-nc-ts"]

        model = GR4JCN(tempfile.mkdtemp())
        
        start_date=dt.datetime(2000,1,1)
        end_date=dt.datetime(2002,1,1)
        
        model.rvi.start_date=start_date
        model.rvi.end_date=dt.datetime(2000,1,2)
        model.rvi.run_name = "test"

        model.rvh.name = "Salmon"
        model.rvh.area = "4250.6"
        model.rvh.elevation = "843.0"
        model.rvh.latitude = 54.4848
        model.rvh.longitude = -123.3659

        model.rvp.params = model.params(0.529, -3.396, 407.29, 1.072, 16.9, 0.947)
        model([ts,])
        rvc = model.outputs["solution"]
        
        for i in range(100):
            xa = assimilateQobsSingleDay(model,rvc,ts,np.array([start_date+dt.timedelta(days=1+i)]),number_members=25,precip_std=0.30,temp_std=2.0,qobs_std=0.15)
            #TODO : NEED TO UPDATE rvc FILE WITH xa