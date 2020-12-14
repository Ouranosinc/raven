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


def assimilateQobsSingleDay(model, xa, ts, days, std, solutions):

    model = deepcopy(model)
    tmp = Path(tempfile.mkdtemp())

    number_members = len(solutions)

    p_fn = tmp / "perturbed_forcing.nc"

    perturbed = {}
    dists = {
        "pr": "gamma",
        "rainfall": "gamma",
        "prsn": "gamma",
    }  # Keyed by standard_name

    model.setup_model_run(
        [
            ts,
        ]
    )  # Force the model to reupdate its ts file. If we don<t do this, the model will use the previous run's data.

    for key, s in std.items():
        nc = model.rvt.get(key)
        with xr.open_dataset(nc.path) as ds:
            da = ds.get(nc.var_name).sel(time=days)
            perturbed[nc.var_name] = perturbation(
                da, dists.get(key, "norm"), s, members=number_members
            )
            if nc.var_name is "qobs":
                qobs_full = da.values
                qobs = da.isel(time=-1).values
    perturbed = xr.Dataset(perturbed)

    # preset the assimilation flag. By default, we will do it.
    do_assimilation = True

    # However, if there are nans (especially in observed flow, but others will make the process fail too), don't assimilate and return the states as-is.
    if perturbed.isnull().any():
        do_assimilation = False  # We will simply run the model with the current rvcs for the given time-step.

    perturbed.to_netcdf(p_fn, mode="w")

    # Generate my list of initial states with replacement of soil0 and soil1

    soil0 = []
    soil1 = []
    qsim_matrix = []
    solutions_end = []

    # Not able to parallelize, TODO later. For each member:
    for i in range(number_members):
        model.rvc.parse(
            solutions[i]
        )  # parse and update model with the member's initial states from previous run
        model.rvc.hru_state = model.rvc.hru_state._replace(
            soil0=xa[0, i], soil1=xa[1, i]
        )  # update the assimilated variables
        model(
            p_fn, nc_index=i, overwrite=True
        )  # run the model with the good timeseries and good ncindex
        soil0.append(
            model.storage.isel(time=-1)["Soil Water[0]"].values
        )  # get variable soil0 at end of run
        soil1.append(
            model.storage.isel(time=-1)["Soil Water[1]"].values
        )  # get variable soil1 at end of run
        qsim_matrix.append(
            model.q_sim.values
        )  # get and stack the members in a Qsim matrix
        solutions_end.append(
            model.outputs["solution"].read_text()
        )  # read the output file for the next period's assimilation

    # stack the variables to be assimilated into a single matrix
    x_matrix = np.column_stack((np.array(soil0), np.array(soil1))).T
    qsim_matrix = np.array(qsim_matrix).squeeze()  # make this an array to display later
    qsim_vector = qsim_matrix[:, -1]  # get last period Qsim for assimilation

    # if there is no problem related to missing Qobs or other vars:
    if do_assimilation:
        # prepare the perturbed Qobs and Qobs errors for each member
        qobs_pert = perturbed["qobs"].isel(time=-1).values
        qobs_error = np.tile(qobs, (1, number_members)) - qobs_pert.T

        # actually do the assimilation, return "xa", the assimilated states for the next period
        xa = update_state(
            x_matrix, qobs_pert.reshape(1, -1), qobs_error, qsim_vector.reshape(1, -1)
        )
    else:

        # If no assimilation, then the unassimilated states are returned.
        xa = x_matrix

    return [xa, qsim_matrix, np.array(qobs_full), model, solutions_end]

    # Set new state variables in the rvc file
    # Updated RVC file after assimilation, ready for next simulation day.


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
    qobs_pert : ndarray (1, n_members)
      Perturbed observed streamflows.
    qobs_error : ndarray (1, n_members)
      Perturbation added to qobs to get qobs_pert.
    qsim : ndarray (1, n_members)
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
    n_members = qobs_pert.shape[1]
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
