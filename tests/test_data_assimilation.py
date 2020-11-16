import datetime as dt
import os
import tempfile
import numpy as np
import xarray as xr
import pytest
import matplotlib.pyplot as plt
from copy import deepcopy
from raven.utilities.data_assimilation import assimilateQobsSingleDay, perturbation
from raven.models import GR4JCN
from .common import TESTDATA
import pdb


def test_perturbation():
    ts = TESTDATA["raven-gr4j-cemaneige-nc-ts"]
    ds = xr.open_dataset(ts)

    tmax = ds.tmax.isel(time=slice(0, 10))
    p_tmax = perturbation(tmax, "norm", .01, members=50)
    np.testing.assert_allclose(p_tmax.mean("members"), tmax, rtol=.1)

    rain = ds.rain.isel(time=slice(30, 60))
    p_rain = perturbation(rain, "gamma", 0.01, members=50)
    np.testing.assert_allclose(p_rain.mean("members"), rain, rtol=.1)

    assert p_tmax.attrs == ds.tmax.attrs
    assert p_rain.attrs == ds.rain.attrs


class TestAssimilationGR4JCN:
    def test_simple(self):
        ts = TESTDATA["raven-gr4j-cemaneige-nc-ts"]
        number_members=25
        std = {"rainfall": 0.30,"prsn": 0.30, "tasmin": 2.0, "tasmax": 2.0, "water_volume_transport_in_river_channel": 0.15}
        model = GR4JCN(tempfile.mkdtemp())

        start_date=dt.datetime(2000,6,1)
        end_date=dt.datetime(2002,6,1)

        model.rvi.start_date=start_date
        model.rvi.end_date=dt.datetime(2000,6,10)
        model.rvi.run_name = "test"

        model.rvh.name = "Salmon"
        model.rvh.area = "4250.6"
        model.rvh.elevation = "843.0"
        model.rvh.latitude = 54.4848
        model.rvh.longitude = -123.3659

        model.rvp.params = model.params(104.1,-1.4127,167.79,4.3798,1.9337,0.56398) #SALMON
        #model.rvp.params = model.params(0.529, -3.396, 407.29, 1.072, 16.9, 0.947)

        model([ts,])
        rvc = model.outputs["solution"]
        model.rvc.parse(rvc.read_text())

        xa = [model.storage["Soil Water[0]"].isel(time=-1).values, model.storage["Soil Water[1]"].isel(time=-1).values]
        xa = np.tile(xa, (number_members, 1)).T

        # Do first run at 7 days
        date_list = [start_date + dt.timedelta(days=x) for x in range(10)]

        [xa,q_assim,q_openloop,model] = assimilateQobsSingleDay(model,xa,ts,date_list,std,number_members=number_members)


        solutions = model.solution
       
        
        for i in range(1,20):
            start_date=start_date+dt.timedelta(days=10)
            date_list = [start_date + dt.timedelta(x) for x in range(10)]
            [xa,q_a,q_ol,model] = assimilateQobsSingleDay(model,xa,ts,date_list,std,number_members=number_members, solutions=solutions)
            model=deepcopy(model)
            rvc=model.outputs["solution"]
            model.rvc.parse(rvc.read_text())
            q_assim=np.concatenate((q_assim,q_a),1)
            q_openloop=np.concatenate((q_openloop,q_ol),1)

        plt.plot(q_assim.transpose(),'r')
        plt.plot(q_openloop.transpose(),'b')
        plt.show()
