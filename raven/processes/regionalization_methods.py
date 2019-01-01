"""
regionalization_method.py

Provide various methods for hydrological regionalization.
"""
try:
    from regionalization_tools import MLR
except ImportError:
    from .regionalization_tools import MLR
import numpy as np
from netCDF4 import Dataset
import datetime as dt
from pywps import Service
from pywps.tests import assert_response_success
from . common import client_for, TESTDATA, CFG_FILE, get_output, urlretrieve

from raven.models import HMETS
import logging
LOGGER = logging.getLogger("PYWPS")
from raven.processes import RavenHMETSProcess, RavenGR4JCemaNeigeProcess, RavenMOHYSEProcess

import pdb


def MLR_parameters_regionalization(ncfilename,start_date, end_date, latitude, longitude, elevation, area,
                                   hydrological_model, inputs_file,
                                   ungauged_properties,
                                   model_parameters,
                                   gauged_properties,
                                   properties_to_use):
    """
    MLR_parameters_regionalization

    INPUTS:
        hydrological_model -- String of the selected hydrological model
        inputs_file -- The ungauged catchment inputs file to work with
        ungauged_properties -- Properties of the ungauged catchment
        model_parameters -- Selected hydrological model calibrated parameters
        gauged_properties -- Properties of all gauged catchments available
        properties_to_use -- The catchments' properties to use

    Function to compute MLR parameters regionalization method.
    """

    MLR_parameters, r_squared = MLR(ungauged_properties,
                                    model_parameters,
                                    gauged_properties,
                                    properties_to_use)

    # Preallocate empty array with final dimensions
    nc_file = Dataset(inputs_file[0], 'r')
    nDays=nc_file.variables["rain"][:].shape[0]
    Qsim = np.zeros(shape=(nDays,1))
    
    # Make parameter vector
    MLR_param=",".join(map(str,MLR_parameters))
        
 
    if hydrological_model == "GR4JCN":
        
        Qs=runLocalGR4JCN(MLR_param,ncfilename,start_date, end_date, latitude, longitude, elevation, area)
        Qsim[0:Qs.shape[0],1]=Qs.flatten() 
        
    elif hydrological_model == "HMETS":
        
        Qs=runLocalHMETS(MLR_param,ncfilename,start_date, end_date, latitude, longitude, elevation, area )
        Qsim[0:Qs.shape[0],1]=Qs.flatten()
            
    elif hydrological_model == "MOHYSE":
        params=",".join(map(str,MLR_parameters[0:7]))
        hrus=",".join(map(str,MLR_parameters[8:9]))
        Qs=runLocalMOHYSE(params,hrus,ncfilename,start_date, end_date, latitude, longitude, elevation, area )
        Qsim[0:Qs.shape[0],1]=Qs.flatten()
        
    else:
        print("The hydrological model selected is not supported for \
              regionalization")

    return Qsim


def distance_only_regionalization(ncfilename,start_date, end_date, latitude, longitude, elevation, area,
                                  hydrological_model,
                                  inputs_file,
                                  number_donors,
                                  model_parameters):
    """
    distance_only_regionalization

    INPUTS:
        hydrological_model -- String of the selected hydrological model
        inputs_file -- The ungauged catchment inputs file to work with
        number_donors -- The number of gauged catchments to use
        model_parameters -- Selected hydrological model calibrated parameters

    Function to compute distance only regionalization method.
    """
    # Preallocate empty array with final dimensions
    nc_file = Dataset(inputs_file[0], 'r')
    nDays=nc_file.variables["rain"][:].shape[0]
    Qsim = np.zeros(shape=(nDays,number_donors))
  #  ts,start_date,end_date,init,name,run_name,area,elevation,latitude,longitude
    
   
    
    for j in range(0, number_donors):
        
        tmp=model_parameters.values[j,2:]
        params=",".join(map(str,tmp))
        
        if hydrological_model == "GR4JCN":
            
            Qs=runLocalGR4JCN(params,ncfilename,start_date, end_date, latitude, longitude, elevation, area )
            Qsim[0:Qs.shape[0],j]=Qs.flatten()       
        
        elif hydrological_model == "HMETS":
            
            Qs=runLocalHMETS(params,ncfilename,start_date, end_date, latitude, longitude, elevation, area )
            Qsim[0:Qs.shape[0],j]=Qs.flatten()
            
        elif hydrological_model == "MOHYSE":
            params=",".join(map(str,tmp[0:7]))
            hrus=",".join(map(str,tmp[8:9]))
            Qs=runLocalMOHYSE(params,hrus,ncfilename,start_date, end_date, latitude, longitude, elevation, area )
            Qsim[0:Qs.shape[0],j]=Qs.flatten()
        
        else:
            print("The hydrological model selected is not supported for regionalization")
 
    Qsim=np.delete(Qsim,list(range(Qs.shape[0],Qsim.shape[0])),axis=0);

    
    return Qsim


def distance_MLR_regionalization(ncfilename,start_date, end_date, latitude, longitude, elevation, area,
                                 hydrological_model, inputs_file,
                                 number_donors,
                                 ungauged_properties,
                                 model_parameters,
                                 gauged_properties,
                                 properties_to_use):
    """
    distance_MLR_regionalization

    INPUTS:
        hydrological_model -- String of the selected hydrological model
        inputs_file -- The ungauged catchment inputs file to work with
        number_donors -- The number of gauged catchments to use
        ungauged_properties -- Properties of the ungauged catchment
        model_parameters -- Selected hydrological model calibrated parameters
        gauged_properties -- Properties of all gauged catchments available
        properties_to_use -- The catchments' properties to use

    Function to compute distance MLR regionalization method.
    """

    MLR_parameters, r_squared = MLR(ungauged_properties,
                                    model_parameters,
                                    gauged_properties,
                                    properties_to_use)

    model_parameters = model_parameters.values[:, 2:]

    for i in range(0, r_squared.shape[0]):
        # if we have an R2 > 0.5 then we consider this to be a better estimator
        if r_squared[i] > 0.5:
            model_parameters[:, i] = MLR_parameters[i]

    nc_file = Dataset(inputs_file[0], 'r')
    nDays=nc_file.variables["rain"][:].shape[0]
    Qsim = np.zeros(shape=(nDays,number_donors))
    

    for j in range(0, number_donors):
        
       params=",".join(map(str,model_parameters[j,:])) 

       
       if hydrological_model == "GR4JCN":
        
           Qs=runLocalGR4JCN(params,ncfilename,start_date, end_date, latitude, longitude, elevation, area)
           Qsim[0:Qs.shape[0],j]=Qs.flatten() 
        
       elif hydrological_model == "HMETS":
        
           Qs=runLocalHMETS(params,ncfilename,start_date, end_date, latitude, longitude, elevation, area )
           Qsim[0:Qs.shape[0],j]=Qs.flatten()
            
       elif hydrological_model == "MOHYSE":
           params=",".join(map(str,model_parameters[j,0:7]))
           hrus=",".join(map(str,model_parameters[j,8:9]))
           
           Qs=runLocalMOHYSE(params,hrus,ncfilename,start_date, end_date, latitude, longitude, elevation, area )
           Qsim[0:Qs.shape[0],j]=Qs.flatten()
        
       else:
           print("The hydrological model selected is not supported for \
                 regionalization")

    return Qsim


def runLocalHMETS(params,ncfilename,start_date, end_date, latitude, longitude, elevation, area ):
    client = client_for(Service(processes=[RavenHMETSProcess(), ], cfgfiles=CFG_FILE))
    
    datainputs = "ts=files@xlink:href=file://{ts};" \
                     "params={params};" \
                     "start_date={start_date};" \
                     "end_date={end_date};" \
                     "init={init};" \
                     "name={name};" \
                     "run_name={run_name};" \
                     "area={area};" \
                     "latitude={latitude};" \
                     "longitude={longitude};" \
                     "elevation={elevation};" \
    .format(ts=ncfilename,
                    params=params,
                    start_date=start_date,
                    end_date=end_date,
                    init='000,000',
                    name='HMETS',
                    run_name='hmets_regionalize',
                    area=area,
                    elevation=elevation,
                    latitude=latitude,
                    longitude=longitude,
                    )
    
    resp = client.get(service='WPS', request='Execute', version='1.0.0', identifier='raven-hmets',datainputs=datainputs)
           
            
    assert_response_success(resp)
    out = get_output(resp.xml)
            
    tmp_file, _ = urlretrieve(out['hydrograph'])
    Qs=Dataset(tmp_file,'r')
    Qs=Qs['q_sim'][:,:].data
    
    
    return Qs

def runLocalGR4JCN(params,ncfilename,start_date, end_date, latitude, longitude, elevation, area ):
    client = client_for(Service(processes=[RavenGR4JCemaNeigeProcess(), ], cfgfiles=CFG_FILE))
    
    datainputs = "ts=files@xlink:href=file://{ts};" \
                     "params={params};" \
                     "start_date={start_date};" \
                     "end_date={end_date};" \
                     "init={init};" \
                     "name={name};" \
                     "run_name={run_name};" \
                     "area={area};" \
                     "latitude={latitude};" \
                     "longitude={longitude};" \
                     "elevation={elevation};" \
    .format(ts=ncfilename,
                    params=params,
                    start_date=start_date,
                    end_date=end_date,
                    init='000,000',
                    name='GR4J',
                    run_name='GR4J_regionalize',
                    area=area,
                    elevation=elevation,
                    latitude=latitude,
                    longitude=longitude,
                    )
    
 
    resp = client.get(service='WPS', request='Execute', version='1.0.0', identifier='raven-gr4j-cemaneige',datainputs=datainputs)
           
            
    assert_response_success(resp)
    out = get_output(resp.xml)
            
    tmp_file, _ = urlretrieve(out['hydrograph'])
    Qs=Dataset(tmp_file,'r')
    Qs=Qs['q_sim'][:,:].data
    
    
    return Qs

def runLocalMOHYSE(params,hrus,ncfilename,start_date, end_date, latitude, longitude, elevation, area ):
    client = client_for(Service(processes=[RavenMOHYSEProcess(), ], cfgfiles=CFG_FILE))
    
    datainputs = "ts=files@xlink:href=file://{ts};" \
                     "params={params};" \
                     "start_date={start_date};" \
                     "end_date={end_date};" \
                     "init={init};" \
                     "name={name};" \
                     "run_name={run_name};" \
                     "area={area};" \
                     "latitude={latitude};" \
                     "longitude={longitude};" \
                     "elevation={elevation};" \
                     "hrus={hrus};"\
    .format(ts=ncfilename,
                    params=params,
                    hrus = hrus,
                    start_date=start_date,
                    end_date=end_date,
                    init='000,000',
                    name='MOHYSE',
                    run_name='MOHYSE_regionalize',
                    area=area,
                    elevation=elevation,
                    latitude=latitude,
                    longitude=longitude,
                    )
    
   
    resp = client.get(service='WPS', request='Execute', version='1.0.0', identifier='raven-mohyse',datainputs=datainputs)
           
            
    assert_response_success(resp)
    out = get_output(resp.xml)
            
    tmp_file, _ = urlretrieve(out['hydrograph'])
    Qs=Dataset(tmp_file,'r')
    Qs=Qs['q_sim'][:,:].data
    
    
    return Qs