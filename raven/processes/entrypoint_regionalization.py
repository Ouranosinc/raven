
"""
entrypoint_regionalization.py

INPUTS:
    hydrological_model -- The hydrological model to use, limited as of now to:
                          'HMETS', 'MOHYSE', 'GR4J'.
    inputs_file -- Ungauged catchment inputs file containing the following fields:
        rain (extracted from boundaries beforehand)
        snow (extracted from boundaries beforehand)
        tasmax (extracted from boundaries beforehand)
        tasmin (extracted from boundaries beforehand)
        surface_area (extracted from GIS beforehand)
        day of year (extracted from dataset used for meteo)
        dates (extracted from dataset used for meteo)
        latitude (centroid, extracted from GIS beforehand)
        longitude (centroid, extracted from GIS beforehand)
        land cover fractions (determined by GIS tools beforehand)
        avg slope (determined by GIS beforehand)
        avg elevation (determined by GIS beforehand)

    regionalization_method -- An integer with the value-argument pairs:
        1- Multiple Linear Regression (MLR)
        2- Spatial Proximity (SP)
        3- Spatial Proximity with Inverse Distance Weighting (IDW) of
           donor basins (SP_IDW)
        4- Spatial Proximity with IDW and regression-based augmentation (RA)
            (SP_IDW_RA)
        5- Physical Similarity (PS)
        6- Physical Similarity with IDW of donor basins (PS_IDW)
        7- Physical Similarity with IDW and RA (PS_IDW_RA)

    number_donors -- A positive integer from [1:10] that allows averaging the
        hydrographs of number_donors donors. For example, a value of 5
        indicates that we take the average hydrograph of the 5 closest or most
        similar donors. No impact when used for MLR. IDW computes the average 
        using inverse weighting.

    min_NSE -- The minimum calibration NSE required for a donor catchment to be
        considered. If the closest or most similar catchment has an NSE below
        the min_NSE, then it is skipped and the next closest/similar catchment
        is selected.

OUTPUT:
    This function returns two arrays. Qsim = [n_days x number_donors] 
    in which the columns are the individual hydrographs used to produce the 
    averaged hydrograph, in order of closest to farthest (or most to least similar)
    The second is "best_estimate", which is the average (or IDW average) of the 
    Qsim hydrographs.
"""

import numpy as np

from regionalization_methods import (
        MLR_parameters_regionalization,
        distance_only_regionalization,
        distance_MLR_regionalization
)
from regionalization_tools import (
        get_array_of_gauged_properties,
        filter_array_by_NSE,
        rank_donor_catchments,
        get_ungauged_properties,
        IDW
)


def Regionalisation(hydrological_model,
                               inputs_file,
                               regionalization_method=2,
                               number_donors=1,
                               min_NSE=0.6):
    
    # Find the properties we wish to use for multiple linear regression (MLR)
    # and similarity.
    properties_to_use = ['centroid_latitude',
                         'centroid_longitude',
                         'surface_area',
                         'avg_elevation',
                         'avg_slope',
                         'land_forest',
                         'land_grass',
                         'land_impervious',
                         'land_urban']

    # Get the ungauged catchment properties from the inputs_file for the
    # regionalization scheme.
    ungauged_properties = get_ungauged_properties(inputs_file,properties_to_use)


    # Load gauged catchments properties and hydrological model parameters.
    gauged_properties, model_parameters = get_array_of_gauged_properties(
                                            hydrological_model)

    # Filter arrays according to the minimum value allowed.
    gauged_properties, model_parameters = filter_array_by_NSE(
                                            gauged_properties,
                                            model_parameters,
                                            min_NSE)

    # Check to see if we have enough data, otherwise raise error
    if (gauged_properties.shape[0] < number_donors and
            regionalization_method >= 2):
        raise ValueError("hydrological_model and minimum NSE threshold \
                         combination is too strict for the number of donor \
                         basins. Please reduce the number of donor basins OR \
                         reduce the minimum NSE threshold")

    # Rank the matrix according to the similarity or distance.
    if (regionalization_method >= 5):  # physical similarity
        regionalization_tag = 'similarity'
    else:
        # Filter by distance if MLR (option 1) is selected, it is faster
        # than by similarity.
        regionalization_tag = 'distance'

    # Create a dataframe that is identical to gauged_properties arrays.
    gauged_properties, model_parameters = rank_donor_catchments(
                                            ungauged_properties,
                                            model_parameters,
                                            gauged_properties,
                                            properties_to_use,
                                            regionalization_tag
                                            )

    if regionalization_method == 1:
        # Multiple linear regression (MLR)
        Qsim = MLR_parameters_regionalization(hydrological_model,
                                              inputs_file,
                                              ungauged_properties,
                                              model_parameters,
                                              gauged_properties,
                                              properties_to_use)
        best_estimate = Qsim

    elif (regionalization_method == 2 or regionalization_method == 5):
        # Simple averaging
        Qsim = distance_only_regionalization(hydrological_model,
                                             inputs_file,
                                             number_donors,
                                             model_parameters)
        best_estimate = np.average(Qsim, 1)

    elif (regionalization_method == 3 or regionalization_method == 6):
        # Inverse distance weighting (IDW)
        Qsim = distance_only_regionalization(hydrological_model,
                                             inputs_file,
                                             number_donors,
                                             model_parameters)
        best_estimate = IDW(Qsim, number_donors, gauged_properties['distance'])

    elif (regionalization_method == 4 or regionalization_method == 7):
        # Inverse distance weighting (IDW)
        Qsim = distance_MLR_regionalization(hydrological_model,
                                            inputs_file,
                                            number_donors,
                                            ungauged_properties,
                                            model_parameters,
                                            gauged_properties,
                                            properties_to_use)
        best_estimate = IDW(Qsim, number_donors, gauged_properties['distance'])

    else:
        print("Unknown regionalization method")

    return best_estimate, Qsim


"""
These are test codes to be removed upon deploy
"""
best_estimate, Qsim = entrypoint_regionalization("GR4JCN", 'inputData3.nc',1,2,0.6)
