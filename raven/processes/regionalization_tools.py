"""
regionalization_tools.py

Provide various tools for hydrological regionalization.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from netCDF4 import Dataset

import logging

LOGGER = logging.getLogger("PYWPS")


def get_array_of_gauged_properties(hydrological_model):
    """
    get_array_or_gauged_properties

    INPUTS:
        hydrological_model -- String of the selected hydrological model

    This function return an array of gauged catchments properties as well as
    model NSE and optimal parameters of the hydrological model.
    """

    # Load general gauged catchment properties and create dataframe df.
    gauged_properties = pd.read_csv(
        Path(__file__).parent.parent.parent / 'tests/testdata/regionalisation_data/gauged_catchment_properties.csv',
        dtype=float)
    types_dict = {'NAME': str, 'STATE': str, 'COUNTRY': str}
    for col, col_type in types_dict.items():
        gauged_properties[col] = gauged_properties[col].astype(col_type)

    gauged_properties = gauged_properties.drop('ID', 1)
    gauged_properties = gauged_properties.drop('NAME', 1)
    gauged_properties = gauged_properties.drop('STATE', 1)
    gauged_properties = gauged_properties.drop('COUNTRY', 1)

    # Get data for the hydrological_model
    filepath = str(Path(
        __file__).parent.parent.parent) + '/tests/testdata/regionalisation_data/' + hydrological_model + '_parameters.csv'
    model_parameters = pd.read_csv(filepath, dtype=float)
    model_parameters = model_parameters.drop('ID', 1)

    # Make sure that the matrices are the same length, else raise error.
    if np.size(gauged_properties, 0) != np.size(model_parameters, 0):
        raise ValueError("gauged_properties and model_parameters for \
                         hydrological model " + hydrological_model +
                         " do not have the same number of gauged catchments. \
                         Check for mismatch.")

    return gauged_properties, model_parameters


def filter_array_by_NSE(gauged_properties, model_parameters, min_NSE):
    """
    filter_array_by_NSE

    INPUTS:
        gauged_properties -- Properties of all gauged catchments available
        model_parameters -- Hydrological model calibrated parameters
        min_NSE -- Minimum Nash-Sutcliffe (NSE) value to select a catchment

    Function to filter the arrays of gauged catchments properties and model
    parameters according to minimum NSE required by the user.
    """

    gauged_properties_filtered = \
        gauged_properties.loc[model_parameters['NASH'] >= min_NSE]
    model_parameters_filtered = \
        model_parameters.loc[model_parameters['NASH'] >= min_NSE]
    gauged_properties_filtered.reset_index(inplace=True)
    model_parameters_filtered.reset_index(inplace=True)

    return gauged_properties_filtered, model_parameters_filtered


def rank_donor_catchments(ungauged_properties,
                          model_parameters,
                          gauged_properties,
                          properties_to_use,
                          rank_type='distance'):
    """
    rank_donor_catchments

    INPUTS:
        model_parameters -- Selected hydrological model calibrated parameters
        ungauged_properties -- Properties of the ungauged catchment
        rank_type -- Regionalization method for the ranking of the catchments
        properties_to_use -- The catchments' properties to use if
                             rank_type='Similarity'

    Function to sort arrays according to the type of sorting (distance or
    similarity) selected by the user.
    """

    # Calculate distance between ungauged and database of gauged catchments.
    dist = np.empty([gauged_properties.shape[0], 1])

    # For the spatial distance:
    if rank_type == 'distance':
        # Calculate the spatial distance between the donor catchments and the
        # ungauged catchment

        for i in range(0, gauged_properties.shape[0]):
            dist[i] = ((gauged_properties.CENTROID_LONGITUDE[i] -
                        ungauged_properties.CENTROID_LONGITUDE) ** 2 +
                       (gauged_properties.CENTROID_LATITUDE[i] -
                        ungauged_properties.CENTROID_LATITUDE) ** 2) ** 0.5

        # Add the physical distance to the dataframe
        gauged_properties = gauged_properties.assign(distance=dist)

    # For the physical similarity:
    elif rank_type == 'similarity':
        # Get the values from the arrays, convert to numpy array
        similar_properties = gauged_properties[properties_to_use].values
        similar_targets = np.tile(
            ungauged_properties[properties_to_use].values,
            (similar_properties.shape[0], 1))

        # Calculate deltas in the database
        delta_max = np.amax(similar_properties, axis=0, keepdims=True)
        delta_min = np.amin(similar_properties, axis=0, keepdims=True)
        deltas = delta_max - delta_min

        # Compute similarity distance
        dist = np.sum(abs(similar_properties - similar_targets) / deltas, 1)
        gauged_properties = gauged_properties.assign(distance=dist)

    else:
        raise ValueError("Regionalization method " + rank_type +
                         " not recognized")

    # The distance is now in the vectors. It has to be sorted, keep the index
    # and also sort the model_parameters matrix
    gauged_properties = gauged_properties.sort_values(by=['distance'])
    model_parameters = model_parameters.reindex(gauged_properties.index)
    gauged_properties.reset_index(drop=True, inplace=True)
    model_parameters.reset_index(drop=True, inplace=True)

    return gauged_properties, model_parameters


def get_ungauged_properties(inputs_file, properties_to_use, latitude, longitude):
    """
    get_ungauged_properties

    INPUTS:
        inputs_file -- The ungauged catchment inputs file to work with
        properties_to_use -- The catchments' properties to use

    Function to get the ungauged catchment's properties from the netCDF file.
    """

    # Read the netCDF file
    nc_file = Dataset(inputs_file[0], 'r')

    """
    # For now, until we test for real and pass the parameter
    properties_to_use = ['CENTROID_LATITUDE',
                         'CENTROID_LONGITUDE',
                         'AREA',
                         'avg_elevation',
                         'avg_slope',
                         'land_forest',
                         'land_grass',
                         'land_impervious',
                         'land_urban']
    """
    # Defile the ungauged catchment properties
    ungauged_properties = pd.DataFrame({'CENTROID_LATITUDE': [latitude], 'CENTROID_LONGITUDE': [longitude]})

    return ungauged_properties


def IDW(Qsim, number_donors, distance):
    """
    IDW - Inverse Distance Weighting

    INPUTS:
        Qsim -- Simulated streamflow
        number_donors -- The number of gauged catchments to use
        distance -- The distance between the ungauged and the gauged catchments

    Function to calculate and weight the Qsim by the distance factor
    """

    # In IDW, weights are 1 / distance
    weights = 1.0 / distance[0:number_donors]

    # Make weights sum to one
    weights /= weights.sum(axis=0)

    # Calculate weighted average
    best_estimate = np.dot(Qsim, weights)

    return best_estimate


def MLR(ungauged_properties,
        model_parameters,
        gauged_properties,
        properties_to_use):
    """
    MLR - Multiple Linear Regression

    INPUTS:
        ungauged_properties -- Properties of the ungauged catchment
        model_parameters -- Selected hydrological model calibrated parameters
        gauged_properties -- Properties of all gauged catchments available
        properties_to_use -- The catchments' properties to use

    Function to compute regionalization using MLR.
    """

    # Define X predictors and Y targets for parameters
    properties_to_use_with_const = properties_to_use[:]
    # Add the constant column for the intercept
    properties_to_use_with_const.insert(0, 'const')
    X = gauged_properties[properties_to_use]
    # Add the constant 1 for the ungauged catchment predictors
    predictors = sm.add_constant(ungauged_properties, prepend=True, has_constant='add')

    # Add constants to the gauged predictors
    X = sm.add_constant(X)
    MLR_parameters = []
    r_squared = []

    # Convert to matrix
    model_parameters = model_parameters.iloc[:, :].values

    # Loop each parameter
    for j in range(2, model_parameters.shape[1]):
        # Define predictands and fit statistical model
        Y = model_parameters[:, j]
        regression_model = sm.OLS(Y, X).fit()

        # Collect R-squared (adjusted) and estimated parameters for the gauged
        # catchment
        r_squared.append(regression_model.rsquared_adj)
        tmp = regression_model.predict(exog=predictors)
        MLR_parameters.append(tmp[0])

    # Convert to arrays
    MLR_parameters = np.array(MLR_parameters)
    r_squared = np.array(r_squared)

    return MLR_parameters, r_squared
