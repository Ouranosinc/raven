"""
Tools for hydrological forecasting
"""

from raven.models import get_model
import logging
import xarray as xr
import numpy as np
import pandas as pd
import datetime as dt
import tempfile

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
    rvc = m.outputs['solution']

    return rvc


# Do the actual forecasting step
def perform_forecasting_step(rvc,model,**kwds):
    '''
    Function that might be useful eventually to do a forecast from a model setup.
    '''
    # kwds includes 'ts', the forecast timeseries data
    # Setup the model
    m = get_model(model)()

    # Force the initial conditions
    m.resume(rvc)

    # Set the parameters, start dates, etc. required to run the model and run
    m(overwrite=True, **kwds)

    return m.q_sim


def perform_climatology_esp(model_name,ts,forecast_date, lead_time, **kwds):
    '''
    This function takes the model setup and name as well as forecast data and lead time and returns
    an ESP forecast netcdf. The data comes from the climatology data and thus there is a mechanism
    to get the correct data from the time series and exclude the current year.

    Parameters:
    -----------
    model_name : {'HMETS', 'MOHYSE', 'GR4J', 'HBVEC'}
        model name to instatiate Raven model
    ts : string
        path to the netcdf timeseries to run the model, which should be a historical timeseries.
    forecast_date : datetime datetime object
        date of the forecast issue.
    lead_time : integer
        number of days of forecast, forward looking.
    kwds : Raven model configuration parameters.

    Returns:
    --------
    qsims: Arraty of streamflow values from the ESP method.

    '''
    # Get the timeseries
    tsnc=xr.open_dataset(ts)

    # Prepare model instance
    m = get_model(model_name)()

    # Now find the periods of time for warm-up and forecast and add to the model keywords as the defaults are failing (nanoseconds datetimes do not like the year 0001...)
    start_date=pd.to_datetime(tsnc['time'][0].values)
    start_date=start_date.to_pydatetime()
    kwds['start_date']=start_date

    # Check to make sure forecast date is not in the first year as we need model warm-up.
    dateLimit = start_date.replace(year=start_date.year + 1)
    if (forecast_date.month==2 and forecast_date.day==29):
        forecast_date.replace(day=28)
    if dateLimit>forecast_date:
        print(' FORECAST DATE IS IN THE WARM-UP WINDOW. PLEASE SELECT ANOTHER FORECAST DATE')

    # Now prepare the met data for the run. It will have the data from day 0 to forecast date at first,
    # Then we replace the leadTime values with a cutout from a previous year, run model and repeat.

    # Base meteo block from begining of time series to beginning of forecast date, to run the model and save the initial conditions.
    baseMet=tsnc.sel(time=slice(start_date,forecast_date-dt.timedelta(days=1)))
    tmp_path=tempfile.mkdtemp() # Eventualy we'll write the ts file here.
    qsims = []

    # list of unique years in the dataset:
    avail_years = list(np.unique(tsnc['time.year'].data))

    # Remove the year that we are forecasting. Or else it's cheating!
    avail_years.remove(forecast_date.year)

    # Get the list of variables in the netcdf dataset
    keylist=list(tsnc.keys())

    # Run model to get rvc file after warm-up using base meteo.
    baseMet.to_netcdf(tmp_path + '/climatologyESP.nc')
    kwds['ts']=tmp_path + '/climatologyESP.nc'
    rvc = get_raven_states(model_name, **kwds)

    # Preassign the variables in the "block" dataset. Cheap way to build a dataset object.
    block_ini=tsnc.sel(time=slice(forecast_date,forecast_date-dt.timedelta(days=1)+dt.timedelta(days=lead_time)))

    # Forcing start and end dates here because the default 0001 year is not working with the datetime64 with nanoseconds, only yeras 1642-2256 or whatever are possible.
    kwds['start_date']=pd.to_datetime(block_ini['time'][0].values).to_pydatetime()
    kwds['end_date']=pd.to_datetime(block_ini['time'][-1].values).to_pydatetime()

    # We will iterate this for all forecast years
    for years in avail_years:

        # Use try statement here as it is possible the final year will fail or some variables are not used.
        try:
            # For each available variable in the dataset
            for k in keylist:

                # Ensure it's a timeseries!
                if 'time' in tsnc[k].coords:

                    # Here we extract the data from the timeseries file for the desired date and place it directly following the forecast date so the timeseries are coherent.
                    block=tsnc[k].sel(time=slice(forecast_date.replace(year=years),forecast_date.replace(year=years) + dt.timedelta(days=lead_time-1))).data
                    block_ini[k].loc[dict(time=slice(forecast_date,forecast_date-dt.timedelta(days=1)+dt.timedelta(days=lead_time)))]=block

            # Now we have finished updating "block_ini" so we can use that as the timeseries to RAVEN.
            block_ini.to_netcdf(tmp_path + '/forecastMember.nc')
            kwds['ts'] = tmp_path + '/forecastMember.nc'


            # This works in Raven build 251+. Run the model and get the simulated hydrograph associated to the forecast block
            m.resume(rvc)
            m(overwrite=True, **kwds)

            # Add member to the ensemble.
            qsims.append(m.q_sim.copy(deep=True))

        except:

            # depending on the forecast length, this might happen! Also if the netcdf file has vars such as qobs that don't cover the whole period, this will stop it.
            print('Last available year not used as dataset does not cover the lead time OR bad variable')

    # Concatenate the members through a new dimension for the members and remove unused dims.
    qsims = xr.concat(qsims, dim='member')
    qsims=qsims.squeeze()
    return qsims
