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

import pdb



def assimilateQobsSingleDay(model,rvc,ts,days,number_members=25,precip_std=0.30,temp_std=2.0,qobs_std=0.15):
    tmp = Path(tempfile.mkdtemp())
    # Extract the Qobs, precip and temperature data for the day of interest
    nc_inputs = xr.open_dataset(ts)
    
    nc_inputs=nc_inputs.sel(time=days)
    #pr=nc_inputs['pr'].sel(time=slice(days[0],days[1]))
    #tasmax=nc_inputs['tasmax'].sel(time=slice(days[0],days[1]))
    #tasmin=nc_inputs['tasmin'].sel(time=slice(days[0],days[1]))
    #qobs=nc_inputs['qobs'].sel(time=slice(days[0],days[1]))
    
    # Check to see if all data is OK for assimilation (i.e. is Qobs=Nan or something like that, then skip)
    do_assimilation=True
    for varname, da in nc_inputs.data_vars.items():
        if np.isnan(da.data):
            do_assimilation=False # We will simply run the model with the current rvcs for the given time-step.
    
    # Extract vars individually  
    qobs=nc_inputs['qobs'].data   
    pr=nc_inputs['rain'].data+nc_inputs['snow'].data
    tasmax=nc_inputs['tmax'].data
    tasmin=nc_inputs['tmin'].data    
    
    # Apply sampling input data uncertainty by perturbation (number_members, precip_std, temp_std, qobs_std)
    pr_pert,tasmax_pert,tasmin_pert,qobs_pert,qobs_error=applyHydrometPerturbations(pr,tasmin,tasmax,qobs,number_members, precip_std, temp_std, qobs_std)
    
    # Run the model number_member times, writing a NetCDF each time with adjusted precip/temp
    x_matrix=[]
    qsim_matrix=[]
    for i in range(0,number_members):
        x, qsim = runModelwithShortMeteo(model,rvc,pr,tasmax,tasmin,days,tmp)
        x_matrix.append(x)
        qsim_matrix.append(qsim)
    
    # Apply the actual EnKF on the states to generate new, better states.
    if do_assimilation:
        xa=applyAssimilationOnStates(x,qobs_pert,qobs_error,qsim_matrix)
    else:
        xa=x
    
    return xa    
    
    # Set new state variables in the rvc file
    # Updated RVC file after assimilation, ready for next simulation day.    
    
def applyHydrometPerturbations(pr,tasmax,tasmin,qobs,number_members, precip_std, temp_std, qobs_std):
    
    
    
    # Temperature: Sample from normal distribution
    random_noise_temperature = np.random.normal(0, temp_std, [1,number_members])
    tasmax_pert=tasmax+random_noise_temperature
    tasmin_pert=tasmin+random_noise_temperature
    
    # Qobs: Sample from normal distribution
    qobs_pert=np.random.normal(loc=np.tile(qobs,(1,number_members)), scale=qobs_std*np.tile(qobs,(1,number_members)))
    qobs_error=np.tile(qobs,(1,number_members))-qobs_pert # might have to copy/cat/repmat qobs here.
    
    # Precipitation: Use a gamma function. Use shape and scale parameters independently
    shape_k=(pr**2)/(precip_std*pr)**2
    scale_t=((precip_std*pr)**2) / pr
    pr_pert=np.random.gamma(shape=np.tile(shape_k,(1,number_members)),scale=np.tile(scale_t,(1,number_members)))
    
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
    number_members=qobs_pert.shape[0]
    z=(qsim)*np.ones((number_members,1))
    HA=(qsim)-(z*np.ones((1,number_members)))/number_members
    Y=qobs_pert-qsim    
    Re=(qobs_error*qobs_error.transpose())/number_members # ensemble covariance calculated by equation 19 in: The Ensemble Kalman Filter: theoretical formulation and practical implementation, Evensen 2003, https://link.springer.com/article/10.1007%2Fs10236-003-0036-9
    P=Re+(HA*HA.transpose())/(number_members-1)
    M=(P**-1)*Y
    Z=(HA.transpose())*M
    A=X-((X*np.ones((number_members,1)))*np.ones((1,number_members)))/number_members
    Xa=X+(A*Z)/(number_members-1)
    Xa=np.maximum(Xa,0)
    
    return Xa

def runModelwithShortMeteo(model,rvc,pr,tasmax,tasmin,days,tmp):

    
    pdb.set_trace()
    
    days=np.array(days)
    pr=np.array(pr)
    
    ds = xr.Dataset({
    'pr': xr.DataArray(
                data   = pr,
                dims   = ['time'],
                coords = {'time': days},
                attrs  = {
                    '_FillValue': -999.9,
                    'units'     : 'mm'
                    }
                ),
    'tmax': xr.DataArray(
                data   = tasmax,   # enter data here
                dims   = ['time'],
                coords = {'time': days},
                attrs  = {
                    '_FillValue': -999.9,
                    'units'     : 'degC'
                    }
                ),
    'tmin': xr.DataArray(
                data   = tasmin,   # enter data here
                dims   = ['time'],
                coords = {'time': days},
                attrs  = {
                    '_FillValue': -999.9,
                    'units'     : 'degC'
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
    
    # update the rvc file
    model.resume(rvc)
    
    # run the model with the input files
    model([tsfile, ])
    
    # allocate space and fill in with the state variables
    x=np.empty(2)
    x[0]=model.storage["Soil Water[0]"]
    x[1]=model.storage["Soil Water[1]"]
    qsim=model.q_sim
    
    return x, qsim

