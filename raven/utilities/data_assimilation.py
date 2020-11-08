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



def assimilateQobsSingleDay(model,xa,ts,days,number_members=25,precip_std=0.30,temp_std=2.0,qobs_std=0.15):
    tmp = Path(tempfile.mkdtemp())
    
    # Extract the Qobs, precip and temperature data for the day of interest
    nc_inputs = xr.open_dataset(ts)
    nc_inputs=nc_inputs.sel(time=days)

    
    # Check to see if all data is OK for assimilation (i.e. is Qobs=Nan or something like that, then skip)
    do_assimilation=True
   
    for varname, da in nc_inputs.data_vars.items():
        if np.isnan(da.data).any():
            do_assimilation=False # We will simply run the model with the current rvcs for the given time-step.
    
    # Extract vars individually  
    qobs=np.array([nc_inputs['qobs'].data])   
    pr=np.array([nc_inputs['rain'].data+nc_inputs['snow'].data])
    tasmax=np.array([nc_inputs['tmax'].data])
    tasmin=np.array([nc_inputs['tmin'].data])  

    # Apply sampling input data uncertainty by perturbation (number_members, precip_std, temp_std, qobs_std)
    pr_pert,tasmax_pert,tasmin_pert,qobs_pert,qobs_error=applyHydrometPerturbations(pr,tasmin,tasmax,qobs[0,-1],number_members, precip_std, temp_std, qobs_std)
    
    # Run the model number_member times, writing a NetCDF each time with adjusted precip/temp
    x_matrix=[]
    qsim_matrix=[]
    rvc_main=model.outputs["solution"]
    
    for i in range(0,number_members):
        model.resume(rvc_main)
        
        # THIS IS WHERE IT ALL FAILS.
        model.rvc.hru_state=model.rvc.hru_state._replace(soil0=xa[0,i])
        model.rvc.hru_state=model.rvc.hru_state._replace(soil1=xa[1,i])
        
        x, qsim = runModelwithShortMeteo(model,pr_pert[:,i],tasmax_pert[:,i],tasmin_pert[:,i],days,tmp)
        x_matrix.append(x)
        qsim_matrix.append(qsim.data)
    
    qsim_matrix=np.array(qsim_matrix)

    qsim_matrix=qsim_matrix[:,:,0]
    
 
    x_matrix=np.array(x_matrix)
    x_matrix=x_matrix[:,:,-1].transpose()
    # Apply the actual EnKF on the states to generate new, better states.
    if do_assimilation:
        xa=applyAssimilationOnStates(x_matrix,qobs_pert,qobs_error,qsim_matrix[:,-1])
    else:
        xa=x_matrix

    return [xa,qsim_matrix,qobs]    
    
    # Set new state variables in the rvc file
    # Updated RVC file after assimilation, ready for next simulation day.    
    
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
    model([tsfile, ])
    
    # allocate space and fill in with the state variables
    x=[]
    x.append(model.storage["Soil Water[0]"])
    x.append(model.storage["Soil Water[1]"])
    qsim=model.q_sim
   
    
    return x, qsim

