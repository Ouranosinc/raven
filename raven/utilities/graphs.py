"""
Library to perform graphs for the streamflow time series analysis.

The following graphs can be plotted:
    - hydrograph
    - mean_annual_hydrograph
    - spaghetti_annual_hydrograph

"""

from matplotlib import pyplot as plt
import xarray as xr
import pandas as pd
import numpy as np


def hydrograph(file_list):
    """
    annual_hydrograph

    INPUTS:
        file_list -- Raven output files containing simulated streamflows

    Create a graphic of the hydrograph for each model simulation.
    """

    ds = [xr.open_dataset(file) for file in file_list]
  
    # Get time data for the plot
    dates = pd.DatetimeIndex(ds[0].time.values)
    first_date = dates.min().strftime('%Y/%m/%d')
    last_date = dates.max().strftime('%Y/%m/%d')

    basin_name = ds[0].basin_name.values[0]  # selected basin name

    fig, ax = plt.subplots()  # initialize figure

    # Plot the observed streamflows if available
    if hasattr(ds[0], 'q_obs'):
        q_obs = ds[0].q_obs
        plt.plot(
            dates,
            q_obs,
            linewidth=2,
            label='obs')

    # Plot the simulated streamflows for each hydrological model
    q_sim = [h.q_sim for h in ds]
    for sim in q_sim:
        plt.plot(
            dates,
            sim,
            linewidth=2,
            label='sim: ' + '<model_name>')

    plt.xlim([first_date, last_date])
    plt.ylim(bottom=0, top=None)
    ax.set_xlabel('Time')
    ax.set_ylabel(r'$Streamflow [m^3s^{{-1}}]$')
    ax.set_title('Hydrograph between ' + first_date + ' and ' + last_date +
                 '\n Selected basin: ' + basin_name + ' basin.')
    ax.legend()
    ax.grid()

    plt.xticks(rotation=90)
    plt.tight_layout()

    return fig


def mean_annual_hydrograph(file_list):
    """
    mean_annual_hydrograph

    INPUTS:
        file_list -- Raven output files containing simulated streamflows

    Create a graphic of the mean hydrological cycle for each model simulation.
    """

    # Time series for the plot
    ds = [xr.open_dataset(file) for file in file_list]

    # Get time data for the plot
    dates = pd.DatetimeIndex(ds[0].time.values)
    first_date = dates.min().strftime('%Y/%m/%d')
    last_date = dates.max().strftime('%Y/%m/%d')
   
    basin_name = ds[0].basin_name.values[0]  # selected basin name

    fig, ax = plt.subplots()  # initialize figure

    # Plot the observed streamflows if available
    if hasattr(ds[0], 'q_obs'):
        q_obs = ds[0].q_obs
        mah_obs = q_obs.groupby('time.dayofyear').mean()

        plt.plot(
            mah_obs.dayofyear,
            mah_obs,
            linewidth=2,
            label='obs')

    # Plot the simulated streamflows for each hydrological model
    q_sim = [h.q_sim for h in ds]
    mahs_sim = [q.groupby('time.dayofyear').mean() for q in q_sim]

    for mah in mahs_sim:
        plt.plot(
            mah.dayofyear,
            mah,
            linewidth=2,
            label='sim: ' + '<model_name>')

    plt.xticks(
        np.linspace(0, 365, 13)[:-1], (
            'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
            'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec',
            ),
        )

    plt.xlim(0, mah_obs.shape[0])
    plt.ylim(bottom=0, top=None)
    ax.set_xlabel('Time')
    ax.set_ylabel(r'$Streamflow [m^3s^{{-1}}]$')
    ax.set_title('Hydrograph between ' + first_date + ' and ' + last_date +
                 '\n Selected basin: ' + basin_name + ' basin.')
    ax.legend()
    ax.grid()

    plt.tight_layout()

    return fig


def spaghetti_annual_hydrograph(file):
    """
    spaghetti_annual_hydrograph

    INPUTS:
        file -- Raven output file containing simulated streamflows of one model

    Create a spaghetti plot of the mean hydrological cycle for one model
    simulations. The mean simulation is also displayed.
    """

    # Time series for the plot
    ds = xr.open_dataset(file)
    
    # Get time data for the plot
    dates = pd.DatetimeIndex(ds.time.values)
    first_date = dates.min().strftime('%Y/%m/%d')
    last_date = dates.max().strftime('%Y/%m/%d')

    basin_name = ds.basin_name.values[0]  # selected basin name

    fig, ax = plt.subplots()  # initialize figure

    # Plot the observed streamflows if available
    if hasattr(ds, 'q_obs'):
        q_obs = ds.q_obs
        mah_obs = q_obs.groupby('time.year')
        mah_obs_mean = q_obs.groupby('time.dayofyear').mean()

        for year in mah_obs.groups.keys():
            plt.plot(
                np.arange(1, q_obs.values[mah_obs.groups[year]].shape[0]+1, 1),
                q_obs.values[mah_obs.groups[year]],
                linewidth=1,
                color='C0')

        plt.plot(
                mah_obs_mean.dayofyear,
                mah_obs_mean,
                linewidth=2,
                color='C0',
                label='obs')

    # Plot the simulated streamflows for each hydrological model
    q_sim = ds.q_sim
    mah_sim = q_sim.groupby('time.year')
    mah_sim_mean = q_sim.groupby('time.dayofyear').mean()

    for year in mah_sim.groups.keys():
        plt.plot(
            np.arange(1, q_sim.values[mah_sim.groups[year]].shape[0]+1, 1),
            q_sim.values[mah_sim.groups[year]],
            linewidth=1,
            color='C1')

    plt.plot(
        mah_sim_mean.dayofyear,
        mah_sim_mean,
        linewidth=2,
        color='C1',
        label='sim: ' + '<model_name>')

    plt.xticks(
        np.linspace(0, 365, 13)[:-1], (
            'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
            'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec',
            ),
        )

    plt.xlim(0, mah_sim_mean.shape[0])
    plt.ylim(bottom=0, top=None)
    ax.set_xlabel('Time')
    ax.set_ylabel(r'$Streamflow [m^3s^{{-1}}]$')
    ax.set_title('Spaghetti annual hydrograph between ' + first_date +
                 ' and ' + last_date + '\n Selected basin: ' + basin_name +
                 ' basin.')
    ax.legend()
    ax.grid()

    plt.tight_layout()

    return fig
