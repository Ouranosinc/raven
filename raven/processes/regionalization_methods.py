"""
regionalization_method.py

Provide various methods for hydrological regionalization.
"""

from raven.processes.regionalization_tools import MLR
import numpy as np
from netCDF4 import Dataset

from raven.models import HMETS, GR4JCemaneige, MOHYSE

import logging

LOGGER = logging.getLogger("PYWPS")


def MLR_parameters_regionalization(ncfilename, start_date, end_date, latitude, longitude, elevation, area,
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
    nDays = nc_file.variables["rain"][:].shape[0]
    Qsim = np.zeros(shape=(nDays, 1))

    # Make parameter vector
    MLR_param = MLR_parameters

    if hydrological_model == "GR4JCN":

        Qs = runLocalGR4JCN(MLR_param, ncfilename, start_date, end_date, latitude, area)
        Qsim[0:Qs.shape[0], 1] = Qs

    elif hydrological_model == "HMETS":

        Qs = run_local_HMETS(MLR_param, ncfilename, start_date, end_date, latitude, area)
        Qsim[0:Qs.shape[0], 1] = Qs

    elif hydrological_model == "MOHYSE":

        Qs = run_local_MOHYSE(MLR_param, ncfilename, start_date, end_date, latitude, area)
        Qsim[0:Qs.shape[0], 1] = Qs

    else:
        print("The hydrological model selected is not supported for \
              regionalization")

    return Qsim


def distance_only_regionalization(ncfilename, start_date, end_date, latitude, longitude, elevation, area,
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
    nDays = nc_file.variables["rain"][:].shape[0]
    Qsim = np.zeros(shape=(nDays, number_donors))

    for j in range(0, number_donors):

        params = model_parameters.values[j, 2:]

        if hydrological_model == "GR4JCN":

            Qs = runLocalGR4JCN(params, ncfilename, start_date, end_date, latitude, area)
            Qsim[0:Qs.shape[0], j] = Qs

        elif hydrological_model == "HMETS":

            Qs = run_local_HMETS(params, ncfilename, start_date, end_date, latitude, area)
            Qsim[0:Qs.shape[0], j] = Qs

        elif hydrological_model == "MOHYSE":

            Qs = run_local_MOHYSE(params, ncfilename, start_date, end_date, latitude, area)
            Qsim[0:Qs.shape[0], j] = Qs

        else:
            print("The hydrological model selected is not supported for regionalization")

    Qsim = np.delete(Qsim, list(range(Qs.shape[0], Qsim.shape[0])), axis=0);

    return Qsim


def distance_MLR_regionalization(ncfilename, start_date, end_date, latitude, longitude, elevation, area,
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
    nDays = nc_file.variables["rain"][:].shape[0]
    Qsim = np.zeros(shape=(nDays, number_donors))

    for j in range(0, number_donors):

        params = model_parameters[j, :]

        if hydrological_model == "GR4JCN":

            Qs = runLocalGR4JCN(params, ncfilename, start_date, end_date, latitude, area)
            Qsim[0:Qs.shape[0], j] = Qs

        elif hydrological_model == "HMETS":

            Qs = run_local_HMETS(params, ncfilename, start_date, end_date, latitude, area)
            Qsim[0:Qs.shape[0], j] = Qs

        elif hydrological_model == "MOHYSE":

            Qs = run_local_MOHYSE(params, ncfilename, start_date, end_date, latitude, area)
            Qsim[0:Qs.shape[0], j] = Qs

        else:
            print("The hydrological model selected is not supported for \
                 regionalization")

    return Qsim


# LOCAL MODEL DEFS
def run_local_HMETS(params, ts, start_date, end_date, latitude, area):
    hmets = HMETS(workdir='/tmp/test')

    params = hmets.RVP.params(GAMMA_SHAPE=params[0],
                              GAMMA_SCALE=params[1],
                              GAMMA_SHAPE2=params[2],
                              GAMMA_SCALE2=params[3],
                              MIN_MELT_FACTOR=params[4],
                              MAX_MELT_FACTOR=params[5],
                              DD_MELT_TEMP=params[6],
                              DD_AGGRADATION=params[7],
                              SNOW_SWI_MIN=params[8],
                              SNOW_SWI_MAX=params[9],
                              SWI_REDUCT_COEFF=params[10],
                              DD_REFREEZE_TEMP=params[11],
                              REFREEZE_FACTOR=params[12],
                              REFREEZE_EXP=params[13],
                              PET_CORRECTION=params[14],
                              HMETS_RUNOFF_COEFF=params[15],
                              PERC_COEFF=params[16],
                              BASEFLOW_COEFF_1=params[17],
                              BASEFLOW_COEFF_2=params[18],
                              TOPSOIL=params[19],
                              PHREATIC=params[20])

    hmets([ts, ], overwrite=True, rvp={'params': params}, rvi={'start_date': start_date, 'end_date': end_date},
          rvh={'area': area, 'latitude': latitude})
    hmets.diagnostics['DIAG_RMSE']
    Qs = hmets.hydrograph['q_sim'].data
    Qs = Qs.flatten()

    return Qs


def runLocalGR4JCN(params, ts, start_date, end_date, latitude, area):
    gr4j = GR4JCemaneige(workdir='/tmp/test')

    params = gr4j.RVP.params(GR4J_X1=params[0],
                             GR4J_X2=params[1],
                             GR4J_X3=params[2],
                             GR4J_X4=params[3],
                             CEMANEIGE_X1=params[4],
                             CEMANEIGE_X2=params[5])

    gr4j([ts, ], overwrite=True, rvp={'params': params}, rvi={'start_date': start_date, 'end_date': end_date},
         rvh={'area': area, 'latitude': latitude})
    gr4j.diagnostics['DIAG_RMSE']
    Qs = gr4j.hydrograph['q_sim'].data
    Qs = Qs.flatten()
    return Qs


def run_local_MOHYSE(params, ts, start_date, end_date, latitude, area):
    parSet = MOHYSE.RVP.params(par_x01=params[0],
                               par_x02=params[1],
                               par_x03=params[2],
                               par_x04=params[3],
                               par_x05=params[4],
                               par_x06=params[5],
                               par_x07=params[6],
                               par_x08=params[7])

    hrus = MOHYSE.RVH.hrus(par_x09=params[8], par_x10=params[9])

    mohyse = MOHYSE(workdir='/tmp/test')
    mohyse([ts, ], overwrite=True, rvp={'params': parSet}, rvi={'start_date': start_date, 'end_date': end_date},
           rvh={'area': area, 'hrus': hrus, 'latitude': latitude})
    mohyse.diagnostics['DIAG_RMSE']
    Qs = mohyse.hydrograph['q_sim'].data
    Qs = Qs.flatten()

    return Qs
