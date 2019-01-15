
"""
regionalization_tools.py

Provide various tools for hydrological regionalization.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
import xarray as xr
from raven.models import get_model

import logging
LOGGER = logging.getLogger("PYWPS")
DATA_DIR = Path(__file__).parent.parent.parent / 'tests' / 'testdata' / 'regionalisation_data'


def regionalize(method, model, latitude, longitude, size=5, min_NSE=0.6, properties=None, **kwds):
    """Perform regionalization for catchment whose outlet is defined by coordinates.

    Parameters
    ----------
    method : {'MLR', 'SP', 'PS', 'SP_IDW', 'PS_IDW', 'SP_IDW_RA', 'PS_IDW_RA'}
      Name of the regionalization method to use.
    model : {'HMETS', 'GR4JCN', 'MOHYSE'}
      Model name.
    latitude : float
      Coordinate of the catchment's outlet.
    longitude : float
      Coordinate of the catchment's outlet.
    size : int
      Number of catchments to use in the regionalization.
    min_NSE : float
      Minimum calibration NSE value required to be considered as a donor.
    properties : sequence
      Name of catchment properties to include in analysis. Defaults to all.

    Returns
    -------
    (qsim, ensemble)


    """
    # TODO: Include list of available properties in docstring.

    # Get the ungauged catchment properties from the inputs_file for the regionalization scheme.
    props = get_ungauged_properties(latitude, longitude)
    kwds.update(props)
    ungauged_properties = pd.DataFrame(props, index=['target'])

    if properties is None:
        properties = list(props.keys())
    else:
        ungauged_properties = ungauged_properties[properties]

    # Load gauged catchments properties and hydrological model parameters.
    gauged_prop = read_gauged_properties().drop(columns=['NAME', 'STATE', 'COUNTRY'])
    gauged_prop = gauged_prop[properties]
    gauged_nash, gauged_params = read_gauged_params(model)

    # Filter on NSE
    valid = gauged_nash > min_NSE
    filtered_params = gauged_params.where(valid).dropna()
    filtered_prop = gauged_prop.where(valid).dropna()

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
    reg_params = regionalization_params(method, sparams, sprop, ungauged_properties)

    # Run the model over all parameters and create ensemble DataArray
    m = get_model(model)
    qsims = []
    for params in reg_params:
        kwds['params'] = params
        m.run(overwrite=True, **kwds)
        qsims.append(m.hydrograph)

    qsims = xr.concat(qsims, 'ens')

    # 3. Aggregate runs into a single result
    if method in ['MLR', 'SP', 'PS']:  # Average (one realization for MLR, so no effect).
        qsim = qsims.mean(dim='ens')
    elif 'IWD' in method:
        qsim = IDW(qsims, dist)
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

    a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2

    c = 2 * np.arcsin(np.sqrt(a))
    km = 6367 * c
    return km


def distance(gauged, ungauged):
    """Return geographic distance between ungauged and database of gauged catchments."""
    lon, lat = ungauged.longitude, ungauged.latitude
    lons, lats = gauged.latitude, gauged.latitude

    return pd.Series(data=haversine(lons.values, lats.values, lon.values, lat.values), index=gauged.index)


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


def regionalization_params(method, gauged_params, gauged_properties, ungauged_properties):
    """Return the model parameters to use for the regionalization.

    Parameters
    ----------
    method : {'MLR', 'SP', 'PS', 'SP_IDW', 'PS_IDW', 'SP_IDW_RA', 'PS_IDW_RA'}
      Name of the regionalization method to use.


    Returns
    -------
    list
      List of model parameters to be used for the regionalization.
    """

    if method == 'MLR' or 'RA' in method:
        mlr_params, r2 = multiple_linear_regression(gauged_properties, gauged_params, ungauged_properties)

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

    # Calculate weighted average
    return qsims.dot(weights, dims='ens')


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

    mlr_parameters = [r.predict(exog=predictors)[0] for r in regression]
    r2 = [r.rsquared_adj for r in regression]

    return mlr_parameters, r2


