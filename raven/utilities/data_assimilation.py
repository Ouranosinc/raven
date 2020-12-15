# -*- coding: utf-8 -*-
"""
Created on Sun Nov  1 20:48:03 2020

@author: Richard
"""


"""
model = Raven model instance, preset with parameters etc.
xa = the set of state variables. In this case, soil0 and soil1 from GR4JCN. Will eventually need to update for other models, other variables
ts = the timeseries of inputs needed by model (tas, pr, qobs, etc.)
days = data for the days of assimilation. Assimilation is performed after last day.
number_members = number of EnKF members
std = variables for uncertainty estimation of input and streamflow variables.
    ex: std = {"rainfall": 0.30,"prsn": 0.30, "tasmin": 2.0, "tasmax": 2.0, "water_volume_transport_in_river_channel": 0.15}
precip = standard deviation used to sample precip, uses gamma distribution (fraction of observed value)
temperature = standard deviation used to sample temperature, normal dist (degrees Celcius)
qobs = standard deviation used to sample observed streamflow, normal dist (fraction of observed value)
"""

import xarray as xr
import numpy as np
import tempfile
from pathlib import Path
from copy import deepcopy
from raven.models.state import HRU_NC_MAP

"""
Suggestion:
This was initially suggested to make parallel and faster. However, I was blocked and could not make it work
I've reverted to coding it in series but it will need to be made parallel as right now it's way too slow (4h for 300 assimilation periods)

Then use model.rvt.nc_index to run parallel simulations with different hru_state.

model(p_fn, hru_state=[list of hru_states to start the simulation with], nc_index=range(n_members))

You can probably then do something like

last_storage = xr.concat(model.storage, dim="members").isel(time=-1)
last_storage["Soil Water[0]"]

"""


def assimilate(model, ts, q_obs, keys, basin_states, hru_states, days):
    """Assimilate streamflow over one day.

    Parameters
    ----------
    model : raven.Model
      Raven model instance configured to run.
    ts : str, Path
      Perturbed time series.
    keys : tuple
      Name of hru_state attributes to be assimilated, for example ("soil0", "soil1").
    basin_states : sequence
      Model initial conditions, BasinStateVariables instances.
    hru_states : sequence
      Model initial conditions, HRUStateVariables instances.
    ts : str, Path, list
      Input netCDF file names.
    days : datetime
      Dates ...
    std : dict
      Standard deviation of the perturbation noise for each input variable, keyed by variable standard name.

    """
    qkey = "water_volume_transport_in_river_channel"

    if len(basin_states) != len(hru_states):
        raise ValueError("`basin_states` and `hru_states` must have the same length.")

    model = deepcopy(model)

    # Number of members
    n_members = len(basin_states)

    # Run simulation with perturbed inputs
    model(ts, hru_state=hru_states, basin_state=basin_states, nc_index=range(n_members))

    # Extract final states (n_states, n_members)
    f_hru_states, f_basin_states = model.get_final_state()
    x_matrix = np.array([[getattr(state, key) for key in keys] for state in f_hru_states]).T

    # Sanity check
    if x_matrix.shape != (len(keys), n_members):
        raise ValueError

    #vnames = [HRU_NC_MAP[k] for k in keys]
    #x_matrix = model.storage.isel(time=-1)[vnames].to_array()

    # Last time step for assimilation
    qsim_vector = model.q_sim.isel(time=-1, nbasins=0).values

    # If there are problems related to missing Qobs or other variables, do not assimilate.
    with xr.open_dataset(ts) as perturbed:
        if perturbed.isnull().any().to_array().any():
            xa = x_matrix
        else:
            # Prepare the perturbed Qobs and Qobs errors for each member
            qobs_pert = perturbed[qkey].isel(time=-1)
            qobs_error = q_obs.isel(time=-1) - qobs_pert

            # Do the assimilation, return the assimilated states for the next period
            xa = update_state(x_matrix, qobs_pert.values, qobs_error.values, qsim_vector)

    return [xa, model]


def perturbation(da, dist, std, **kwargs):
    """Return perturbed time series.

    Parameters
    ----------
    da : DataArray
      Input time series to be perturbed.
    dist : {"norm", "gamma"}
      Name of statistical distribution from which random perturbations are drawn.
    std : float
      Standard deviation of random perturbation.
    kwargs : dict
      Name and size of additional dimensions, apart from time.

    """
    nt = len(da.time)
    size = list(kwargs.values()) + [nt]
    dims = list(kwargs.keys()) + ["time"]

    if dist == "norm":
        r = np.random.normal(0, std, size=size)
        out = da + xr.DataArray(r, dims=dims, coords={"time": da.time})

    elif dist == "gamma":
        shape = (da ** 2) / (std * da) ** 2
        scale = ((std * da) ** 2) / da
        r = np.nan_to_num(np.random.gamma(shape=shape, scale=scale, size=size), nan=0.0)
        out = xr.DataArray(r, dims=dims, coords={"time": da.time})

    else:
        raise AttributeError(f"{dist} is not supported.")

    out.attrs.update(da.attrs)
    return out


def update_state(x, qobs_pert, qobs_error, qsim):
    """
    Update model state by assimilation.

    Parameters
    ----------
    x : ndarray (n_states, n_members)
      Model state initial values.
    qobs_pert : ndarray (n_members)
      Perturbed observed streamflows.
    qobs_error : ndarray (n_members)
      Perturbation added to qobs to get qobs_pert.
    qsim : ndarray (n_members)
      Simulated streamflows.

    Returns
    -------
    ndarray (n_states, n_members)
      Model state values after assimilation.

    Reference
    ---------
    The Ensemble Kalman Filter: theoretical formulation and practical implementation, Evensen 2003
    https://link.springer.com/article/10.1007%2Fs10236-003-0036-9

    University of colorado report on the efficient implementation of EnKF, Jan Mandel, 2006
    http://ccm.ucdenver.edu/reports/rep231.pdf
    """
    n_states, n_members = np.shape(x)

    # Make sure arrays have shape (1, n_members)
    qobs_pert = np.atleast_2d(qobs_pert)
    qobs_error = np.atleast_2d(qobs_error)
    qsim = np.atleast_2d(qsim)

    z = np.dot(qsim, np.ones((n_members, 1)))
    ha = qsim - (z * np.ones((1, n_members))) / n_members
    y = qobs_pert - qsim

    # Equations 4.1 from Mandel, 2006
    re = np.dot(qobs_error, qobs_error.transpose()) / n_members
    p = re + (np.dot(ha, ha.transpose())) / (n_members - 1)
    m = np.dot(p ** -1, y)
    z = (ha.transpose()) * m
    a = x - np.dot((np.dot(x, np.ones((n_members, 1)))), np.ones((1, n_members))) / n_members
    xa = x + (np.dot(a, z)) / (n_members - 1)
    xa = np.maximum(xa, 0)

    return xa
