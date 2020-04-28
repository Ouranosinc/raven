import pytest
import datetime as dt
import numpy as np
import xarray as xr
import json
import netCDF4 as nc
from pywps import Service
from pywps.tests import assert_response_success
import os

from . common import client_for, TESTDATA, CFG_FILE, get_output, urlretrieve
from raven.processes import ClimatologyEspProcess
import pdb



class TestClimatologyESP:
    def test_simple(self):
        client = client_for(Service(processes=[ClimatologyEspProcess(), ], cfgfiles=CFG_FILE))
        
        params = '9.5019, 0.2774, 6.3942, 0.6884, 1.2875, 5.4134, 2.3641, 0.0973, 0.0464, 0.1998, 0.0222, -1.0919, ' \
                 '2.6851, 0.3740, 1.0000, 0.4739, 0.0114, 0.0243, 0.0069, 310.7211, 916.1947'         

        #Date of the forecast that will be used to determine the members of the climatology-based ESP (same day of year of all other years)
        forecast_date=dt.datetime(1956,7,1)
        lead_time=60 # Number of days for lead time
        model = 'HMETS'
        
        
        datainputs = "ts=files@xlink:href=file://{ts};" \
                     "params={params};" \
                     "init={init};" \
                     "name={name};" \
                     "run_name={run_name};" \
                     "area={area};" \
                     "latitude={latitude};" \
                     "longitude={longitude};" \
                     "elevation={elevation};" \
                     "forecast_date={forecast_date};" \
                     "lead_time={lead_time};" \
                     "model_name={model_name};" \
            .format(ts=TESTDATA['climatologyESP'],
                    params=params,
                    init='155,455',
                    name='Salmon',
                    run_name='test-climatologyESP',
                    area='4250.6',
                    elevation='843.0',
                    latitude=54.4848,
                    longitude=-123.3659,
                    forecast_date=forecast_date,
                    lead_time=lead_time,
                    model_name=model,
                    )
        
        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='climatology_esp',
            datainputs=datainputs)
        pdb.set_trace()
        assert_response_success(resp)
        out = get_output(resp.xml)
        assert 'ensemble' in out
