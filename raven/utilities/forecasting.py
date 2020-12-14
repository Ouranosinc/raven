#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 09:11:58 2020

@author: ets
"""

"""
Tools for hydrological forecasting
"""

import datetime as dt
import logging
import re
import tempfile
import warnings
from pathlib import Path

import fiona
import numpy as np
import pandas as pd
import fiona
from zipfile import ZipFile

import rioxarray
import xarray as xr
from xclim import subset
from raven.models import get_model

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
        msg = (
            "Forecast date is within the warm-up period. Select another forecast date."
        )
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
            msg = (
                f"Year {years} has been removed because it is the last year in the dataset and does not cover the "
                f"forecast duration."
            )
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


def get_hindcast_day(region_coll,date, climate_model='GEPS'):
    '''
    This function generates a forecast dataset that can be used to run raven. 
    Data comes from the CASPAR archive and must be aggregated such that each file
    contains forecast data for a single day, but for all forecast timesteps and
    all members.
    
    The code takes the region shapefile, the forecast date required, and the 
    climate_model to use, here GEPS by default, but eventually could be 
    GEPS, GDPS, REPS or RDPS.
    '''

    # Get the file locations and filenames as a function of the climate model and date
    [ds,times]=get_CASPAR_dataset(climate_model,date)
   

    return get_subsetted_forecast(region_coll,ds,times,True)

def get_CASPAR_dataset(climate_model,date):

    if climate_model == 'GEPS':
        file_url='https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/birdhouse/caspar/daily/GEPS_' + dt.datetime.strftime(date,'%Y%m%d') + '.nc'
        ds = xr.open_dataset(file_url)       
        # Here we also extract the times at 6-hour intervals as Raven must have
        # constant timesteps and GEPS goes to 6 hours
        start = pd.to_datetime(ds.time[0].values)
        times = [start + dt.timedelta(hours=n) for n in range(0, 384, 6)]
    else:
        # Eventually: GDPS, RDPS and REPS
        raise NotImplementedError("Only the GEPS model is currently supported")

    # Checking that these exist.
    for f in ["pr", "tas"]:
        if f not in ds:
            raise AttributeError(f"'{f}' not present in dataset")

    return ds, times
        
def get_ECCC_dataset(climate_model):

    if climate_model == "GEPS":
        # Eventually the file will find a permanent home, until then let's use the test folder.
        #file_url = "https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/datasets/forecasts/eccc_geps/GEPS_latest.ncml"
        file_url = "https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/birdhouse/testdata/geps_forecast/GEPS_latest.ncml"
        ds = xr.open_dataset(file_url)
        # Here we also extract the times at 6-hour intervals as Raven must have
        # constant timesteps and GEPS goes to 6 hours
        start = pd.to_datetime(ds.time[0].values)
        times = [start + dt.timedelta(hours=n) for n in range(0, 384, 6)]
    else:
        # Eventually: GDPS, RDPS and REPS
        raise NotImplementedError("Only the GEPS model is currently supported")


    # Checking that these exist. IF the files are still processing, possible that one or both are not available!
    for f in ["pr", "tas"]:
        if f not in ds:
            raise AttributeError(f"'{f}' not present in dataset")

    return ds, times


def get_recent_ECCC_forecast(region_coll, climate_model="GEPS"):
    """
    This function generates a forecast dataset that can be used to run raven.
    Data comes from the ECCC datamart and collected daily. It is aggregated
    such that each file contains forecast data for a single day, but for all
    forecast timesteps and all members.

    The code takes the region shapefile and the climate_model to use, here GEPS
    by default, but eventually could be GEPS, GDPS, REPS or RDPS.

    Parameters
    ----------
    region_coll : fiona.collection.Collection
      The region vectors.
    climate_model : str
      Type of climate model, for now only "GEPS" is supported.

    Returns
    -------
    forecast : xararray.Dataset
      The forecast dataset.

    """

    [ds,times] = get_ECCC_dataset(climate_model)
    
    return get_subsetted_forecast(region_coll,ds,times,False)


def get_subsetted_forecast(region_coll,ds,times,is_caspar):
    '''
    This function takes a dataset, a region and the time sampling array and returns
    the subsetted values for the given region and times
    
    Parameters
    ----------
    region_coll : fiona.collection.Collection
      The region vectors.
    ds : xarray.Dataset
      The dataset containing the raw, worldwide forecast data
    times: dt.datetime
      The array of times required to do the forecast.
    is_caspar: boolean
      True if the data comes from Caspar, false otherwise. Used to define 
      lat/lon on rotated grid.
      
    Returns
    -------
    forecast : xararray.Dataset
      The forecast dataset.
      
    '''
    # Extract the bounding box to subset the entire forecast grid to something
    # more manageable
    lon_min = region_coll.bounds[0]
    lon_max = region_coll.bounds[2]
    lat_min = region_coll.bounds[1]
    lat_max = region_coll.bounds[3]

    # Subset the data to the desired location (bounding box) and times
    ds = subset.subset_bbox(
        ds, lon_bnds=[lon_min, lon_max], lat_bnds=[lat_min, lat_max]
    ).sel(time=times)

    # Rioxarray requires CRS definitions for variables
    # Get CRS, e.g. 4326
    crs = int(re.match("epsg:(\d+)", region_coll.crs["init"]).group(1))

    # Here the name of the variable could differ based on the Caspar file processing
    tas = ds.tas.rio.write_crs(crs)
    pr = ds.pr.rio.write_crs(crs)
    ds = xr.merge([tas, pr])

    # Now apply the mask of the basin contour and average the values to get a single time series
    if is_caspar == True:
        ds.rio.set_spatial_dims('rlon','rlat')   
        ds["rlon"] = ds["rlon"] - 360
        # clip the netcdf and average across space.
        shdf = [region_coll.next()["geometry"]]
        forecast = ds.rio.clip(shdf, crs=crs)
        forecast = forecast.mean(dim={"rlat", "rlon"}, keep_attrs=True)

    else:
        ds.rio.set_spatial_dims("lon", "lat")
        ds["lon"] = ds["lon"] - 360
        # clip the netcdf and average across space.
        shdf = [region_coll.next()["geometry"]]
        forecast = ds.rio.clip(shdf, crs=crs)
        forecast = forecast.mean(dim={"lat", "lon"}, keep_attrs=True)
    



    return forecast

