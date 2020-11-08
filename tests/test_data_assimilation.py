import datetime as dt
import os
import tempfile
import numpy as np
import pytest

from raven.utilities.data_assimilation import assimilateQobsSingleDay
from raven.models import GR4JCN
from .common import TESTDATA
import pdb

class TestAssimilationGR4JCN:
    def test_simple(self):
        ts = TESTDATA["raven-gr4j-cemaneige-nc-ts"]
        number_members=25
        
        model = GR4JCN(tempfile.mkdtemp())
        
        start_date=dt.datetime(2000,1,1)
        end_date=dt.datetime(2002,1,1)
        
        model.rvi.start_date=start_date
        model.rvi.end_date=dt.datetime(2000,1,7)
        model.rvi.run_name = "test"

        model.rvh.name = "Salmon"
        model.rvh.area = "4250.6"
        model.rvh.elevation = "843.0"
        model.rvh.latitude = 54.4848
        model.rvh.longitude = -123.3659

        model.rvp.params = model.params(0.529, -3.396, 407.29, 1.072, 16.9, 0.947)
        model([ts,])
        rvc = model.outputs["solution"]
        
        
        model.resume(rvc)
        xa=[]
        xa.append(model.storage["Soil Water[0]"])
        xa.append(model.storage["Soil Water[1]"])
        xa=np.array(xa)
        xa=np.tile(xa[:,-1],(number_members,1)).transpose()
        
        # Do first run at 7 days
        date_list = [start_date + dt.timedelta(days=x) for x in range(7)]
      
        [xa,q_assim,q_openloop] = assimilateQobsSingleDay(model,xa,ts,date_list,number_members=number_members,precip_std=0.30,temp_std=2.0,qobs_std=0.15)            
        
        for i in range(1,25):
            
            # TODO: Fix this once assimilation works
            date_list = [start_date + dt.timedelta(days=x+7*i) for x in range(7)]
            [xa,q_a,q_ol] = assimilateQobsSingleDay(model,xa,ts,date_list,number_members=number_members,precip_std=0.30,temp_std=2.0,qobs_std=0.15)            
            
            q_assim=np.concatenate((q_assim,q_a),1)
            q_openloop=np.concatenate((q_openloop,q_ol))
        
        print(q_assim)