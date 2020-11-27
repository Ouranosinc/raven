import numpy as np
import pytest
from pywps import Service
from pywps.tests import assert_response_success

import xarray as xr
from raven.processes import ForecastFloodRiskProcess
from .common import client_for, TESTDATA, CFG_FILE, get_output, urlretrieve


class TestForecastEvaluationProcess:

    def test_forecast_floodrisk_deterministic(self):
        client = client_for(Service(processes=[ForecastFloodRiskProcess(), ], cfgfiles=CFG_FILE))
        
        
        datainputs = "fcst=files@xlink:href=file://{fcst};" \
                     "name=fcst;" \
                     "flood_level={flood_level};" \
                     .format(fcst=TESTDATA['floodrisk_det'],flood_level=0.5)
  
        resp = client.get(service='WPS', request='Execute', version='1.0.0', identifier='forecast-floodrisk',datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)
        floodrisk, _ = urlretrieve(out['flood_risk'])
        ds=xr.open_dataset(floodrisk)
        np.testing.assert_almost_equal(ds['fcst'], [1, 0, 0], 2)
        
    def test_forecast_floodrisk_ensemble(self):
        client = client_for(Service(processes=[ForecastFloodRiskProcess(), ], cfgfiles=CFG_FILE))
        
        
        datainputs = "fcst=files@xlink:href=file://{fcst};" \
                     "name=fcst;" \
                     "flood_level={flood_level};" \
                     .format(fcst=TESTDATA['floodrisk_ens'],flood_level=0.5)
  
        resp = client.get(service='WPS', request='Execute', version='1.0.0', identifier='forecast-floodrisk',datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)
        floodrisk, _ = urlretrieve(out['flood_risk'])
        ds=xr.open_dataset(floodrisk)
        np.testing.assert_almost_equal(ds['fcst'], [0.2, 0.2, 0.8], 2)
        
     