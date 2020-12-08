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
        
        # get timeseries
        ts = TESTDATA["raven-gr4j-cemaneige-nc-ts"]
        
        # set number of members. Using 7 here to make it easier to find and debug.
        number_members=7

        # set parameters for the assimilation
        std = {"rainfall": 0.30,"prsn": 0.30, "tasmin": 2.0, "tasmax": 2.0, "water_volume_transport_in_river_channel": 0.15}
        
        # make a GR4JCN model instance
        model = GR4JCN(tempfile.mkdtemp())
        
        # set the start and end dates for the first assimilation period, warm-up
        start_date=dt.datetime(2000,6,1)
        end_date=dt.datetime(2000,6,10)

        # update the rvi / rvh info. Not using kwds here as I want to also force the state vars, etc.
        model.rvi.start_date=start_date
        model.rvi.end_date=end_date
        model.rvi.run_name = "test"

        model.rvh.name = "Salmon"
        model.rvh.area = "4250.6"
        model.rvh.elevation = "843.0"
        model.rvh.latitude = 54.4848
        model.rvh.longitude = -123.3659

        model.rvp.params = model.params(0.1353389,-.005067198, 576.8007, 6.986121, 1.102917, 0.9224778) #SALMON

        # run the model
        model([ts,])
       
        # prebuild an ensemble of "number_members" initial states, so a copy, of the output solutions that are read through read_text. I can update them later with the parse method.
        solutions_end=[]
        for i in range(number_members):
            solutions_end.append(model.outputs['solution'].read_text())
            
        # get and update the desired parameters at the end of the run.
        xa = [model.storage["Soil Water[0]"].isel(time=-1).values, model.storage["Soil Water[1]"].isel(time=-1).values]

        # also need to tile them to get the same number of members as the "solutions".
        xa = np.tile(xa, (number_members, 1)).T

        # Do first run at 10 days to give some time to stabilize and get a bit of differences in the final state variables. generate date list.
        date_list = [start_date + dt.timedelta(days=x) for x in range(10)]

        # perform the assimilation step here
        [xa,q_assim,q_obs,model, solutions_end] = assimilateQobsSingleDay(model,xa,ts,date_list,std,solutions=solutions_end)
       
        #solutions = model.solution #This should be more efficient but I simply cannot make it work. once in memory, can't reuse it and change values within to perform assimilation.
       
        # Update the start-time becasue we just ran the model for 10 days, let's connect to that
        start_date=start_date+dt.timedelta(days=10)


        number_iterations=4 # let's just do 4 iterations for this test case to keep things short
        assim_timestep=3 #3 days between each assimilation is reasonable and a good balance between performance and speed.

        # This is the loop that will perform the assimilation over the time window that we need.        
        for i in range(number_iterations):
            
            # generate a new series of dates to run, from the end of the last period to the end of the new period.
            date_list = [start_date + dt.timedelta(x) for x in range(assim_timestep)]
            
            # update this in the rvi
            model.rvi.start_date=date_list[0]
            model.rvi.end_date=date_list[-1]
            
            # Perform the assimilation. Note that we use the solutions_end from the previous run and return an updated one after assimilation.
            [xa,q_a,q_o,model,solutions_end] = assimilateQobsSingleDay(model,xa,ts,date_list,std, solutions=solutions_end)
            
            # reupdate the new start_date
            start_date=start_date+dt.timedelta(days=assim_timestep)
            
            # concatenate the assimilated and observed flows for this short period to the lot. At the end, this is what will be plotted to see the relative peformance of the assimilation method.            
            q_assim=np.concatenate((q_assim,q_a),1)
            q_obs=np.concatenate((q_obs,q_o),0)
           
        # Now do a run on the whole period to compare the assimilation with an open_loop run (no assimilation, straight simulation)
        model = GR4JCN(tempfile.mkdtemp())

        start_date=dt.datetime(2000,6,1)
        end_date=start_date+dt.timedelta(days=10+number_iterations*assim_timestep-1)
        
        # We've reset the model to a new instance so let's repopulate. Probably not optimal...
        model.rvi.start_date=start_date
        model.rvi.end_date=end_date
        model.rvi.run_name = "test"

        model.rvh.name = "Salmon"
        model.rvh.area = "4250.6"
        model.rvh.elevation = "843.0"
        model.rvh.latitude = 54.4848
        model.rvh.longitude = -123.3659

        model.rvp.params = model.params(0.1353389,-.005067198, 576.8007, 6.986121, 1.102917, 0.9224778) #SALMON

        model([ts,])
        
        # We can now plot everything!
        plt.plot(q_assim.T,'r') # plot the assimilated flows
        plt.plot(q_obs.T,'b') # plot the observed flows
        plt.plot(model.q_sim,'g') # plot the open_loop (simulation with no assimilation)
        plt.show()
        
        assert q_assim.shape[0]==number_members
        assert q_assim.shape[1]==10+number_iterations*assim_timestep