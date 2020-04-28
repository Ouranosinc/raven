from pathlib import Path
import logging

from pywps import ComplexInput, ComplexOutput, LiteralInput
from pywps import FORMATS
from pywps import Format
from pywps import Process
from . import wpsio as wio

import pandas as pd
import xarray as xr
import numpy as np
from raven.models import get_model
import datetime as dt
import tempfile


class ClimatologyEspProcess(Process):
    def __init__(self):
        
        
        fdate=LiteralInput('forecast_date', 'date of the forcast',
                               abstract='Date that the climatology-based ESP ensemble will be performed',
                               data_type='dateTime',
                               min_occurs=1,
                               max_occurs=1)
        leadtime=LiteralInput('lead_time', 'Forecast lead-time',
                               abstract='duration of the forecast in days',
                               data_type='integer',
                               default=30,
                               min_occurs=0,
                               max_occurs=1)
        
        # The number of parameters will depend on the selected model in model_name
        params=LiteralInput('params', 'Comma separated list of model parameters',
                               abstract='Parameters to run the model',
                               data_type='string',
                               min_occurs=1,
                               max_occurs=1)
                  
        inputs = [fdate, leadtime, params, wio.ts, wio.latitude, wio.longitude, wio.name,
              wio.model_name, wio.area, wio.elevation]

        outputs = [wio.forecast]
    
        super(ClimatologyEspProcess, self).__init__(
            self._handler,
            identifier="climatology_esp",
            title="",
            version="1.0",
            abstract="",
            metadata=[],
            inputs=inputs,
            outputs=outputs,
            keywords=[],
            status_supported=True,
            store_supported=True)

    def _handler(self, request, response):

        # Collect info from call
        kwds = {}
        for key, val in request.inputs.items():
            kwds[key] = request.inputs[key][0].data
        
        # Get timeseries data that will be used to do the climatological ESP.
        ts=request.inputs['ts'][0].file
        tsnc=xr.open_dataset(ts)
        
        # Get info from kwds but remove the ones that are not understood by RAVEN.
        # Have to delete these because or else I get an error: no config key named model_name
        forecast_date = kwds['forecast_date']
        del kwds['forecast_date']
        lead_time = kwds['lead_time']
        del kwds['lead_time']
        model_name = kwds['model_name']
        del kwds['model_name']  

        # Get the model parameters, transform them to a list of floats and write them back to the kwds config.
        params=kwds['params']
        csv = params.replace('(', '').replace(')', '')
        params= list(map(float, csv.split(',')))
        kwds['params']=params  
          
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
        
        # Base meteo block from begining of time series to end of forecast date. With a refactor using initial conditions to initiate forecast, this will have to be reduced
        # to the beginning of the time series to the day before the forecsat begins and save the initial conditions.
        baseMet=tsnc.sel(time=slice(start_date,forecast_date-dt.timedelta(days=1)+dt.timedelta(days=lead_time)))     
        tmp_path=tempfile.mkdtemp() # Eventualy we'll write the ts file here.
        qsims = []
                
        # list of unique years in the dataset:
        avail_years = list(np.unique(tsnc['time.year'].data))
        
        # Remove the year that we are forecasting. Or else it's cheating!
        avail_years.remove(forecast_date.year)
        
        # Get the list of variables in the netcdf dataset
        keylist=list(tsnc.keys())
       
        # We wil iterate this for all forecast years
        for years in avail_years:

            # Use try statement here as it is possible the final year will fail
            try:
                # For each available variable in the dataset
                for k in keylist:
                    
                    # Ensure it's a timeseries!
                    if 'time' in tsnc[k].coords:
                        
                        # Going to eventually need to refactor to include the use of initial conditions
                        # 1- Run the model from day 1 to forecast-minus-1
                        # 2- save initial conditions
                        # 3- run individual forecasts from initial conditions
                        
                        block=tsnc[k].sel(time=slice(forecast_date.replace(year=years),forecast_date.replace(year=years) + dt.timedelta(days=lead_time-1))).data   
                        baseMet[k].loc[dict(time=slice(forecast_date,forecast_date-dt.timedelta(days=1)+dt.timedelta(days=lead_time)))]=block
                        
                # Now we have finished updating "baseMet" so we can use that as the timeseries to RAVEN.
                baseMet.to_netcdf(tmp_path + '/climatologyESP.nc')
                kwds['ts'] = tmp_path + '/climatologyESP.nc'
                
                # Forcing start and end dates here because the default 0001 year is not working with the datetime64 with nanoseconds, only yeras 1642-2256 or whatever are possible.
                kwds['end_date']=pd.to_datetime(baseMet['time'][-1].values).to_pydatetime()

                # This works in Raven build 251+. Run the model and get the simulated hydrograph associated to the forecast block   
                m(overwrite=True, **kwds)
                
                # Add member to the ensemble.
                qsims.append(m.q_sim.copy(deep=True))

            except:
                
                # depending on the forecast length, this might happen! Also if the netcdf file has vars such as qobs that don't cover the whole period, this will stop it.
                print('Last available year not used as dataset does not cover the lead time OR bad variable')
                
        # COncatenate the members through a new dimension for the members and remove unused dims.
        qsims = xr.concat(qsims, dim='member') 
        qsims=qsims.squeeze()
        
        # Prepare the forecast netcdf result file and send the path to the results output.
        forecastfile=(tmp_path + '/forecast.nc')
        qsims.to_netcdf(forecastfile)
        response.outputs['forecast'].file = forecastfile
        
        return response
