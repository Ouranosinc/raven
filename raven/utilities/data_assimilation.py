# -*- coding: utf-8 -*-
"""
Created on Sun Nov  1 20:48:03 2020

@author: Richard
"""


'''
model = Raven model instance, preset with parameters etc.
ts = the timeseries needed by model
rvc = Raven RVC file that contains the initial model states (needs to be the same size as number_members)
day = data for the days of assimilation. Assimilation is performed after last day.
number_members = number of EnKF members
precip_std = standard deviation used to sample precip (fraction of observed value)
temp_std = standard deviation used to sample temperature (degrees Celcius)
qobs_std = standard deviation used to sample observed streamflow (fraction of observed value)
'''

import xarray as xr
import math
import numpy as np
import tempfile
import datetime as dt
from pathlib import Path
from copy import deepcopy
import os
from raven.models.rv import RavenNcData

import pdb

"""
Suggestion:


Then use model.rvt.nc_index to run parallel simulations with different hru_state.

model(p_fn, hru_state=[list of hru_states to start the simulation with], nc_index=range(n_members))

You can probably then do something like

last_storage = xr.concat(model.storage, dim="members").isel(time=-1)
last_storage["Soil Water[0]"]

"""


def assimilateQobsSingleDay(model, xa,ts,days,std,number_members=25):

    model=deepcopy(model)

    tmp = Path(tempfile.mkdtemp())

    p_fn = tmp / "perturbed_forcing.nc"

    perturbed = {}
    dists = {"pr": "gamma", "rainfall": "gamma", "prsn": "gamma"}  # Keyed by standard_name

    # Magic and/or dark arts are used to do things here:
    for key, s in std.items():
        nc = model.rvt.get(key)
        with xr.open_dataset(nc.path) as ds:
            da = ds.get(nc.var_name).sel(time=days)
            perturbed[nc.var_name] = perturbation(da, dists.get(key, "norm"), s, members=number_members)
            if nc.var_name is "qobs":
                qobs=da.isel(time=-1).values
    perturbed = xr.Dataset(perturbed)
    
    do_assimilation=True
    for varname, da in perturbed.data_vars.items():
        if np.isnan(da.data).any():
            do_assimilation=False # We will simply run the model with the current rvcs for the given time-step.

    perturbed.to_netcdf(p_fn)

    # Generate my list of initial states with replacement of soil0 and soil1
    inistates=[]
    for i in range(number_members):
        inistates.append(model.rvc.hru_state._replace(soil0=xa[0,i], soil1=xa[1,i]))

    model(p_fn, hru_state=inistates, nc_index=range(number_members))
    
    # This does not work, error with 'time' dimension not being available... but it is there.
    # also, the model.storage only returns one value. I think we might have to run the 25 members
    # independently and save them each one by one... unless there is another trick?
    # the 25 q_sims do work though, I am able to get them directly.
    last_storage = model.storage.isel(time=-1)
    soil0=np.array(last_storage["Soil Water[0]"].values)
    soil1=np.array(last_storage["Soil Water[1]"].values)
    x_matrix=np.column_stack((soil0,soil1)).T
    qsim_matrix=model.q_sim.isel(time=-1).values
    # Here I now need to redefine the inputs to ApplyassimilationOnStates, TODO
    # Apply the actual EnKF on the states to generate new, better states.
    
    if do_assimilation:
        qobs_pert=perturbed['qobs'].isel(time=-1).values
        qobs_error=np.tile(qobs,(1,number_members))-qobs_pert.T
        xa=applyAssimilationOnStates(x_matrix,qobs_pert.reshape(1, -1),qobs_error,qsim_matrix.T)
    else:
        xa=x_matrix

    return [xa,qsim_matrix,qobs,model]

    # Set new state variables in the rvc file
    # Updated RVC file after assimilation, ready for next simulation day.


def perturbation(da, dist, std, **kwargs):
    """Return perturbed time series."""
    nt = len(da.time)
    size = list(kwargs.values()) + [nt, ]
    dims = list(kwargs.keys()) + ["time"]

    if dist == "norm":
        r = np.random.normal(0, std, size=size)
        out = da + xr.DataArray(r, dims=dims, coords={"time": da.time})

    elif dist == "gamma":
        shape = (da ** 2) / (std * da) ** 2
        scale = ((std * da) ** 2) / da
        r = np.nan_to_num(np.random.gamma(shape=shape, scale=scale, size=size), nan=0.0)
        out = xr.DataArray(r, dims=dims, coords={"time": da.time})

    out.attrs.update(da.attrs)
    return out


# Richard see perturbation function above
def applyHydrometPerturbations(pr,tasmax,tasmin,qobs,number_members, precip_std, temp_std, qobs_std):

    # Temperature: Sample from normal distribution
    random_noise_temperature = np.random.normal(0, temp_std, (number_members,tasmax.shape[1]))
    tasmax_pert=np.tile(tasmax,(number_members,1))+random_noise_temperature
    tasmin_pert=np.tile(tasmin,(number_members,1))+random_noise_temperature
    tasmax_pert=tasmax_pert.transpose()
    tasmin_pert=tasmin_pert.transpose()

    # Qobs: Sample from normal distribution
    qobs_pert=np.random.normal(loc=np.tile(qobs,(1,number_members)), scale=qobs_std*np.tile(qobs,(1,number_members)))
    qobs_error=np.tile(qobs,(1,number_members))-qobs_pert # might have to copy/cat/repmat qobs here.

    # Precipitation: Use a gamma function. Use shape and scale parameters independently
    shape_k=(pr**2)/(precip_std*pr)**2
    scale_t=((precip_std*pr)**2) / pr
    shape_k=shape_k.transpose()
    scale_t=scale_t.transpose()
    pr_pert=np.random.gamma(shape=np.tile(shape_k,(1,number_members)),scale=np.tile(scale_t,(1,number_members)))
    pr_pert=np.nan_to_num(pr_pert,nan=0.0)
    # return the data. All values should now have size nDays x number_members
    return pr_pert,tasmax_pert,tasmin_pert,qobs_pert,qobs_error


def applyAssimilationOnStates(X,qobs_pert,qobs_error,qsim):
    '''
    X = states, matrix of nStates x number_members
    qobs_pert = perturbed observed streamflows of size 1 x number_members
    qobs_error = amount of error added to qobs to get qobs_pert , size 1 x number_members
    qsim = vector of simulated streamflows, size 1 x number_members

    Taken from http://ccm.ucdenver.edu/reports/rep231.pdf
    University of colorado report on the efficient implementation of EnKF
    Jan Mandel, 2006
    Series of equations 4.1 is implemented here
    '''
    number_members=qobs_pert.shape[1]
    z=np.dot(qsim,np.ones((number_members,1)))
    HA=(qsim)-(z*np.ones((1,number_members)))/number_members
    Y=qobs_pert-qsim
    Re=np.dot(qobs_error,qobs_error.transpose())/number_members # ensemble covariance calculated by equation 19 in: The Ensemble Kalman Filter: theoretical formulation and practical implementation, Evensen 2003, https://link.springer.com/article/10.1007%2Fs10236-003-0036-9
    P=Re+(np.dot(HA,HA.transpose()))/(number_members-1)
    M=np.dot(P**-1,Y)
    Z=(HA.transpose())*M
    A=X-np.dot((np.dot(X,np.ones((number_members,1)))),np.ones((1,number_members)))/number_members
    Xa=X+(np.dot(A,Z))/(number_members-1)
    Xa=np.maximum(Xa,0)

    return Xa

def runModelwithShortMeteo(model,pr,tasmax,tasmin,days,tmp):
    model=deepcopy(model)
    #pdb.set_trace()

    ds = xr.Dataset({
    'pr': xr.DataArray(
                data   = pr,
                dims   = ['time'],
                coords = {'time': days},
                attrs  = {
                    '_FillValue': -999.9,
                    'units'     : 'mm/d'
                    }
                ),
    'tmax': xr.DataArray(
                data   = tasmax,   # enter data here
                dims   = ['time'],
                coords = {'time': days},
                attrs  = {
                    '_FillValue': -999.9,
                    'units'     : 'deg_C'
                    }
                ),
    'tmin': xr.DataArray(
                data   = tasmin,   # enter data here
                dims   = ['time'],
                coords = {'time': days},
                attrs  = {
                    '_FillValue': -999.9,
                    'units'     : 'deg_C'
                    }
                )
            },
    )

    # write the netcdf file
    tsfile=tmp / 'tmp_assim_met.nc'
    ds.to_netcdf(tsfile)

    # Update the model for this run
    model.rvi.start_date=days[0]
    model.rvi.end_date=days[-1]


    # run the model with the input files
    model(tsfile,overwrite=True)
    # DO TESTS HERE FOR MODEL UPDATE PROBLEM
   
    # allocate space and fill in with the state variables
    x=[]
    x.append(model.storage["Soil Water[0]"].copy(deep=True))
    x.append(model.storage["Soil Water[1]"].copy(deep=True))
    qsim=model.q_sim.copy(deep=True)

    return x, qsim

