"""
Tools for hydrological regionalization
"""

from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
import xarray as xr
from raven.models import get_model

import logging
LOGGER = logging.getLogger("PYWPS")


# Added directory for test data (smaller database wth only 10 donor catchments)
DATA_DIR = Path(__file__).parent.parent.parent / 'tests' / 'testdata' / 'regionalisation_data'/ 'tests'

def regionalize(method, model, nash, params=None, props=None, target_props=None, size=5,
                min_NSE=0.6, **kwds):
    """Perform regionalization for catchment whose outlet is defined by coordinates.

    Parameters
    ----------
    method : {'MLR', 'SP', 'PS', 'SP_IDW', 'PS_IDW', 'SP_IDW_RA', 'PS_IDW_RA'}
      Name of the regionalization method to use.
    model : {'HMETS', 'GR4JCN', 'MOHYSE'}
      Model name.
    nash : pd.Series
      NSE values for the parameters of gauged catchments.
    params : pd.DataFrame
      Model parameters of gauged catchments. Needed for all but MRL method.
    props : pd.DataFrame
      Properties of gauged catchments to be analyzed for the regionalization. Needed for MLR and RA methods.
    target_props : pd.Series or dict
      Properties of ungauged catchment. Needed for MLR and RA methods.
    size : int
      Number of catchments to use in the regionalization.
    min_NSE : float
      Minimum calibration NSE value required to be considered as a donor.
    kwds : {}
      Model configuration parameters, including the forcing files (ts).

    Returns
    -------
    (qsim, ensemble)
    qsim = multi-donor averaged predicted streamflow
    ensemble = ensemble of members based on number of donors.

    """
    # TODO: Include list of available properties in docstring.
    # TODO: Add error checking for source, target stuff wrt method chosen.

    # Select properties based on those available in the ungauged properties DataFrame.
    if isinstance(target_props, dict):
        ungauged_properties = pd.Series(target_props)
    elif isinstance(target_props, pd.Series):
        ungauged_properties = target_props
    elif isinstance(target_props, pd.DataFrame):
        ungauged_properties = target_props.to_series()
    else:
        raise ValueError

    
    # Filter on NSE
    valid = nash > min_NSE
    filtered_params = params.where(valid).dropna()
    filtered_prop = props.where(valid).dropna()
    
    
    # Check to see if we have enough data, otherwise raise error
    if len(filtered_prop) < size and method != 'MLR':
        raise ValueError("Hydrological_model and minimum NSE threshold \
                         combination is too strict for the number of donor \
                         basins. Please reduce the number of donor basins OR \
                         reduce the minimum NSE threshold.")

    # Rank the matrix according to the similarity or distance.
    if method in ['PS', 'PS_IDW', 'PS_IDW_RA']:  # Physical similarity
        dist = similarity(filtered_prop, ungauged_properties)
    else:  # Geographical distance.
        dist = distance(filtered_prop, ungauged_properties)

    # Series of distances for the first `size` best donors
    sdist = dist.sort_values().iloc[:size]

    # Pick the donors' model parameters and catchment properties
    sparams = filtered_params.loc[sdist.index]
    sprop = filtered_prop.loc[sdist.index]
    
    
    # Get the list of parameters to run
    reg_params = regionalization_params(method, sparams, sprop, ungauged_properties,filtered_params,filtered_prop)
  
    # Run the model over all parameters and create ensemble DataArray
    m = get_model(model)
    qsims = []
  
    for params in reg_params:
    
        # GENERATES N TIMES SAME HYDROGRAPH BUT USES DIFFERENT PARAMETERS???, TRIED TO FIX BUT AM UNABLE
        kwds['params'] = params
        m.run(overwrite=True, **kwds)
        qsims.append(m.hydrograph)
    
    qsims = xr.concat(qsims, 'ens')

    # 3. Aggregate runs into a single result
    if method in ['MLR', 'SP', 'PS']:  # Average (one realization for MLR, so no effect).
        qsim = qsims.mean(dim='ens')
    elif 'IDW' in method:
        qsim = IDW(qsims, sdist)
    else:
        raise ValueError('No matching algorithm for {}'.format(method))

    
    return qsim, qsims


def read_gauged_properties():
    """Return table of gauged catchments properties over North America.

    Returns
    -------
    pd.DataFrame
      Catchment properties keyed by catchment ID.
    """
    return pd.read_csv(DATA_DIR / 'gauged_catchment_properties.csv',
                       index_col='ID')


def read_gauged_params(model):
    """Return table of NASH-Stucliffe Efficiency values and model parameters for North American catchments.

    Returns
    -------
    pd.DataFrame
      Nash-Sutcliffe Efficiency keyed by catchment ID.
    pd.DataFrame
      Model parameters keyed by catchment ID.
    """

    params = pd.read_csv(DATA_DIR / '{}_parameters.csv'.format(model),
                         index_col='ID')

    return params['NASH'], params.iloc[:, 1:]


def haversine(lon1, lat1, lon2, lat2):
    """
    Return the great circle distance between two points on the earth.

    Parameters
    ----------
    lon1, lat1 : ndarray
        Longitude and latitude coordinates in decimal degrees.
    lon2, lat2 : ndarray
        Longitude and latitude coordinates in decimal degrees.

    Returns
    -------
    ndarray
      Distance between points 1 and 2 [km].

    """
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat1) * np.cos(lat2) * (np.sin(dlon / 2.0) ** 2)

    c = 2 * np.arcsin(np.sqrt(a))
    km = 6367 * c
    return km


def distance(gauged, ungauged):
    """Return geographic distance [km] between ungauged and database of gauged catchments.

    Parameters
    ----------
    gauged : pd.DataFrame
      Table containing columns for longitude and latitude of catchment's centroid.
    ungauged : pd.Series
      Coordinates of the ungauged catchment.

    """
    lon, lat = ungauged.longitude, ungauged.latitude
    lons, lats = gauged.longitude, gauged.latitude
    

    return pd.Series(data=haversine(lons.values, lats.values, lon, lat), index=gauged.index)


def similarity(gauged, ungauged, kind='ptp'):
    """Return similarity measure between gauged and ungauged catchments.

    Parameters
    ----------
    gauged : DataFrame
      Gauged catchment properties.
    ungauged : DataFrame
      Ungauged catchment properties
    kind : {'ptp', 'std', 'iqr'}
      Normalization method: peak to peak (maximum - minimum), standard deviation, interquartile range.

    """

    stats = gauged.describe()

    if kind == 'ptp':
        spread = stats.loc['max'] - stats.loc['min']
    elif kind == 'std':
        spread = stats.loc['std']
    elif kind == 'iqr':
        spread = stats.loc['75%'] - stats.loc['25%']

    d = ungauged.values - gauged.values
    n = np.abs(d) / spread.values
    return pd.Series(data=n.sum(axis=1), index=gauged.index)


def regionalization_params(method, gauged_params, gauged_properties, ungauged_properties, filtered_params, filtered_prop):
    """Return the model parameters to use for the regionalization.

    Parameters
    ----------
    method : {'MLR', 'SP', 'PS', 'SP_IDW', 'PS_IDW', 'SP_IDW_RA', 'PS_IDW_RA'}
      Name of the regionalization method to use.
    gauged_params
      DataFrame of parameters for donor catchments (size = number of donors)
    gauged_properties
      DataFrame of properties of the donor catchments  (size = number of donors)
    ungauged_properties
      DataFrame of properties of the ungauged catchment (size = 1)
    filtered_params --> Required for MLR
      DataFrame of parameters of all filtered catchments (size = all catchments with NSE > min_NSE)
    filtered_prop --> Required for MLR
      DataFrame of properties of all filtered catchments (size = all catchments with NSE > min_NSE)
    
    Returns
    -------
    list
      List of model parameters to be used for the regionalization.
    """

    if method == 'MLR' or 'RA' in method:
        mlr_params, r2 = multiple_linear_regression(filtered_prop, filtered_params, ungauged_properties.to_frame().T)
        
        if method == 'MLR':  # Return the multiple linear regression parameters.
            out = [mlr_params, ]

        elif 'RA' in method:
            gp = gauged_params.copy()

            for p, r, col in zip(mlr_params, r2, gauged_params):
                # If we have an R2 > 0.5 then we consider this to be a better estimator
             
                if r > 0.5:
                    gp[col] = p
              
            out = gp.values
            
    else:
        out = gauged_params.values
   
  
    return out


def get_ungauged_properties(latitude, longitude):
    """
    Return the properties of the watershed whose outlet location is given.

    Parameters
    ----------
    latitude : float
      Coordinate of the catchment's outlet.
    longitude : float
      Coordinate of the catchment's outlet.

    Returns
    -------
    dict
      Catchment properties: area, mean elevation, centroid latitude and longitude, average slope, ...
    """

    # Read the netCDF file
    # nc_file = Dataset(inputs_file[0], 'r')

    """
    # For now, until we test for real and pass the parameter
    properties_to_use = ['latitude',
                         'longitude',
                         'area',
                         'avg_elevation',
                         'avg_slope',
                         'land_forest',
                         'land_grass',
                         'land_impervious',
                         'land_urban']
    """

    return {'latitude': latitude, 'longitude': longitude}


def IDW(qsims, dist):
    """
    Inverse distance weighting.

    Parameters
    ----------
    qsims : DataArray
      Ensemble of hydrogram stacked along the `ens` dimension.
    dist : pd.Series
      Distance from catchment which generated each hydrogram to target catchment.

    Returns
    -------
    DataArray
      Inverse distance weighted average of ensemble.
    """

    # In IDW, weights are 1 / distance
    weights = xr.DataArray(1.0 / dist, dims='ens')

    # Make weights sum to one
    weights /= weights.sum(axis=0)

    # Calculate weighted average. 
    # I GET ERROR WITH INITIAL CODE:
    # AttributeError: 'Dataset' object has no attribute 'dot'
    # return qsims.dot(weights, dims='ens')
    
    # This seems to work but then the output is not a DataArray anymore... 
    # IDW and non_IDW return different structures with this.
    return qsims.q_sim.dot(weights.reset_index('ens'),dims='ens')
    
def multiple_linear_regression(source, params, target):
    """
    Multiple Linear Regression for model parameters over catchment properties.

    Uses known catchment properties and model parameters to estimate model parameter over an
    ungauged catchment using its properties.

    Parameters
    ----------
    source : DataFrame
      Properties of gauged catchments.
    params : DataFrame
      Model parameters of gauged catchments.
    target : DataFrame
      Properties of the ungauged catchment.


    Returns
    -------
    (mrl_params, r2)
      A named tuple of the estimated model parameters and the R2 of the linear regression.
    """
    # Add constants to the gauged predictors
    x = sm.add_constant(source)
    
    # Add the constant 1 for the ungauged catchment predictors
    predictors = sm.add_constant(target, prepend=True, has_constant='add')

    # Perform regression for each parameter
    regression = [sm.OLS(params[param].values, x).fit() for param in params]
    
    # Perform prediction on each parameter based on the predictors    
    mlr_parameters = [r.predict(exog=predictors)[0] for r in regression]
    
    # Extract the adjusted r_squared value for each parameter
    r2 = [r.rsquared_adj for r in regression]
        
    return mlr_parameters, r2
