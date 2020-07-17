"""
Library to perform graphs for the streamflow time series analysis.

The following graphs can be plotted:
    - hydrograph
    - mean_annual_hydrograph
    - spaghetti_annual_hydrograph

"""

import numpy as np
import pandas as pd
import xarray as xr
from matplotlib import pyplot as plt
from scipy import stats

from raven.utilities.mk_test import mk_test_calc
from xclim.core.units import units2pint


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
            label='sim: ' + basin_name)

    # plt.xlim([first_date, last_date])
    plt.ylim(bottom=0, top=None)
    ax.set_xlabel('Time')
    ax.set_ylabel(r'$Streamflow [m^3s^{{-1}}]$')
    ax.set_title('Hydrograph between {} and {}\nSelected basin: {} basin.'.format(first_date, last_date, basin_name))
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
            label='sim: ' + basin_name)

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
    ax.set_title('Hydrograph between {} and {}\nSelected basin: {} basin.'.format(first_date, last_date, basin_name))
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
                np.arange(1, q_obs.values[mah_obs.groups[year]].shape[0] + 1, 1),
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
            np.arange(1, q_sim.values[mah_sim.groups[year]].shape[0] + 1, 1),
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
    ax.set_title('Spaghetti annual hydrograph between {} and {}'
                 '\n Selected basin: {} basin.'.format(first_date, last_date, basin_name))
    ax.legend()
    ax.grid()

    plt.tight_layout()

    return fig


def ts_graphs(file, trend=True, alpha=0.05):
    """
    graphs for time series statistics

    INPUTS:
        file -- xclim file containing streamflow statistics for one run

    Create a figure with the statistics so one can see a trend in the data
    """

    # Time series for the plot
    ds = xr.open_dataset(file)
    titlename = ds['ts_stats'].description
    values = ds['ts_stats'].values[:]

    values[0] = 2
    values[1] = 1
    values[2] = 3
    values[3] = 3
    values[4] = 3
    values[5] = 5
    values[6] = 6
    values[7] = 7
    values[8] = 10
    values[9] = 4

    # Get time data for the plot
    dates = pd.DatetimeIndex(ds.time.values)
    # first_date = dates.min()#.strftime('%Y/%m/%d')
    # last_date = dates.max()#.strftime('%Y/%m/%d')
    res = None

    if trend:
        res = stats.theilslopes(values, dates)
        trd, h, p, z = mk_test_calc(values, alpha=alpha)
        titlename = titlename + ", Mann-Kendall h=" + str(h) + ", p-value=" + str(np.round(p, 4))

    fig, ax = plt.subplots()
    ax.plot(dates, values, label='time-series index')

    # plt.xlim([first_date, last_date])
    ax.set_xlabel('Time')
    ax.set_ylabel(r'$Streamflow [m^3s^{{-1}}]$')

    ax.set_title(titlename)

    ax.grid()
    plt.tight_layout()

    if trend:
        # TODO: This does not work yet, trying to compute the y-value of a datetime x-axis value * a slope...
        ax.plot(dates, res[1] + res[0] * dates,
                linestyle='--',
                linewidth=2,
                label='Sen\'s Slope')

    ax.legend()

    return fig


def ts_fit_graph(ts, params):
    """Create graphic showing an histogram of the data and the distribution fitted to it.

    Parameters
    ----------
    ts : str
      Path to netCDF file storing the time series.
    params : str
      Path to netCDF file storing the distribution parameters.

    Returns
    -------
    fig
    """
    from xclim.indices.generic import get_dist

    n = ts.nbasins.size
    dist = params.attrs['scipy_dist']

    fig, axes = plt.subplots(n, figsize=(10, 6), squeeze=False)

    for i in range(n):
        ax = axes.flat[i]
        ax2 = plt.twinx(ax)
        p = params.isel(nbasins=i)

        # Plot histogram of time series as density then as a normal count.
        density, bins, patches = ax.hist(ts.isel(nbasins=i).dropna(dim='time'),
                                         alpha=.5, density=True, bins='auto',
                                         label="__nolabel__")
        ax2.hist(ts.isel(nbasins=i).dropna(dim='time'), bins=bins, facecolor=(1, 1, 1, 0.01), edgecolor='gray',
                 linewidth=1, )

        # Plot pdf of distribution
        dc = get_dist(dist)(*params.isel(nbasins=i))
        mn = dc.ppf(.01)
        mx = dc.ppf(.99)
        q = np.linspace(mn, mx, 200)
        pdf = dc.pdf(q)
        ps = ', '.join(["{:.1f}".format(x) for x in p.values])
        ax.plot(q, pdf, '-', label="{}({})".format(params.attrs['scipy_dist'], ps))

        # Labels
        ax.set_xlabel("{} (${:~P}$)".format(ts.long_name, units2pint(ts.units)))
        ax.set_ylabel("Probability density")
        ax2.set_ylabel("Histogram count")

        ax.legend(frameon=False)

    plt.tight_layout()
    return fig
