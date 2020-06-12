import pytest
import datetime as dt
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from pywps import Service
from pywps.tests import assert_response_success


from . common import client_for, TESTDATA, CFG_FILE, get_output, urlretrieve
from raven.processes import ClimatologyEspProcess


class TestClimatologyESP:
    def test_simple(self):
        client = client_for(Service(processes=[ClimatologyEspProcess(), ], cfgfiles=CFG_FILE))
        
        #
        #model = 'HMETS'
        #params = '9.5019, 0.2774, 6.3942, 0.6884, 1.2875, 5.4134, 2.3641, 0.0973, 0.0464, 0.1998, 0.0222, -1.0919, ' \
        #         '2.6851, 0.3740, 1.0000, 0.4739, 0.0114, 0.0243, 0.0069, 310.7211, 916.1947'         

        model = 'GR4JCN'
        params = '0.529, -3.396, 407.29, 1.072, 16.9, 0.947'
         
        #Date of the forecast that will be used to determine the members of the climatology-based ESP (same day of year of all other years)
        forecast_date=dt.datetime(1964,3,1)
        lead_time=365 # Number of days for lead time
        
        
        
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
        
        assert_response_success(resp)
        out = get_output(resp.xml)
        assert 'forecast' in out
        
        
        # Display forecast to show it works

        forecast, _ = urlretrieve(out['forecast'])
        tmp = xr.open_dataset(forecast)
        qfcst=tmp['q_sim'][:].data.transpose()
        plt.plot(qfcst)
        plt.show()


