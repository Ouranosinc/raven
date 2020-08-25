#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 09:11:58 2020

@author: ets
"""

"""
Tools for hydrological forecasting
"""

from raven.models import get_model
import logging
import xarray as xr
import numpy as np
import pandas as pd
import datetime as dt
import warnings

LOGGER = logging.getLogger("PYWPS")


# This function gets model states after running the model (i.e. states at the end of the run).
def get_raven_states(model, **kwds):
    """Get the RAVEN states file (.rvc file) after a model run.

    Parameters
    ----------
    model : {'HMETS', 'GR4JCN', 'MOHYSE', 'HBVEC'}
      Model name.
    kwds : {}
      Model configuration parameters, including the forcing files (ts).

    Returns
    -------
    rvc : {}
      Raven model forcing file

    """
    # Run the model and get the rvc file for future hotstart.
    m = get_model(model)()
    m(overwrite=True, **kwds)
    rvc = m.outputs["solution"]

    return rvc


# Do the actual forecasting step
def perform_forecasting_step(rvc, model, **kwds):
    """
    Function that might be useful eventually to do a forecast from a model setup.
    """
    # kwds includes 'ts', the forecast timeseries data
    # Setup the model
    m = get_model(model)()

    # Force the initial conditions
    m.resume(rvc)

    # Set the parameters, start dates, etc. required to run the model and run
    m(overwrite=True, **kwds)

    return m.q_sim


def perform_climatology_esp(model_name, forecast_date, forecast_duration, **kwds):
    """
    This function takes the model setup and name as well as forecast data and duration and returns
    an ESP forecast netcdf. The data comes from the climatology data and thus there is a mechanism
    to get the correct data from the time series and exclude the current year.

    Parameters
    ----------
    model_name : {'HMETS', 'MOHYSE', 'GR4JCN', 'HBVEC'}
      Model name to instantiate Raven model.
    forecast_date : datetime.datetime
      Date of the forecast issue.
    forecast_duration : int
      Number of days of forecast, forward looking.
    kwds : dict
      Raven model configuration parameters.

    Returns
    -------
    qsims
      Array of streamflow values from the ESP method along with list of member years

    """
    # Get the timeseries
    tsnc = xr.open_dataset(kwds["ts"])

    # Prepare model instance
    m = get_model(model_name)()

    # Now find the periods of time for warm-up and forecast and add to the model keywords as the defaults are failing
    # (nanoseconds datetimes do not like the year 0001...)
    start_date = pd.to_datetime(tsnc["time"][0].values)
    start_date = start_date.to_pydatetime()

    kwds["start_date"] = start_date

    # Forecasting from Feb 29th is not ideal, we will replace with Feb 28th.
    # Should not change much in a climatological forecast.
    if forecast_date.month == 2 and forecast_date.day == 29:
        forecast_date.replace(day=28)

    # Check to make sure forecast date is not in the first year as we need model warm-up.
    # We cannot use timedelta because if the dataset happens to start on a leap
    # year, then the timedelta=365 days will not be robust. (and we cannot use timedelta(years=1)...)
    dateLimit = start_date.replace(year=start_date.year + 1)
    if dateLimit > forecast_date:
        msg = "Forecast date is whithin the warm-up period. Select another forecast date."
        warnings.warn(msg)

    # initialize the array of forecast variables
    qsims = []

    # list of unique years in the dataset:
    avail_years = list(np.unique(tsnc["time.year"].data))

    # Take a copy of the forecast initial date before overwriting in the forecast step.
    forecast_date_main = forecast_date

    # Remove the year that we are forecasting. Or else it's cheating!
    avail_years.remove(forecast_date.year)

    # Update the forecast end-date, which will be the day prior to the forecast date.
    # So forecasts warm-up will be from day 1 in the dataset to the forecast date.
    kwds["end_date"] = forecast_date - dt.timedelta(days=1)

    # Run model to get rvc file after warm-up using base meteo.
    rvc = get_raven_states(model_name, **kwds)

    # We need to check which years are long enough (ex: wrapping years, 365-day forecast starting in
    # September 2015 will need data up to August 2016 at least)
    for years in avail_years:
        if forecast_date.replace(year=years) + dt.timedelta(
            days=forecast_duration - 1
        ) > pd.to_datetime(tsnc["time"][-1].values):
            avail_years.remove(years)
            msg = f"Year {years} has been removed because it is the last year in the dataset and does not cover the " \
                  f"forecast duration."
            warnings.warn(msg)

    # We will iterate this for all forecast years
    for years in avail_years:

        # Replace the forecast period start and end dates with the climatological ESP dates for the
        # current member (year)
        forecast_date = forecast_date.replace(year=years)
        kwds["start_date"] = forecast_date
        kwds["end_date"] = forecast_date + dt.timedelta(days=forecast_duration - 1)

        # Setup the initial states from the warm-up and run the model.
        # Note that info on start/end dates and timeseries are in the kwds.
        m.resume(rvc)
        m(overwrite=True, **kwds)

        # Add member to the ensemble and retag the dates to the real forecast dates
        # (or else we will get dates from the climate dataset that cover all years)
        new_member = m.q_sim.copy(deep=True)
        new_member["time"] = pd.date_range(
            forecast_date_main, periods=forecast_duration
        )
        qsims.append(new_member)

    # Concatenate the members through a new dimension for the members and remove unused dims.
    qsims = xr.concat(qsims, dim="member")
    qsims = qsims.squeeze()

    # Add the number of the forecast year as member ID
    qsims["member"] = (["member"], avail_years)

    return qsims
