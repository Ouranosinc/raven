import os

import numpy as np
import pandas as pd
import xarray as xr
from pywps.tests import WpsClient, WpsTestResponse

TESTS_HOME = os.path.abspath(os.path.dirname(__file__))
CFG_FILE = os.path.join(TESTS_HOME, 'test.cfg')

TESTDATA = {'gr4j-cemaneige':
            {'pr': '{0}'.format(os.path.join(TESTS_HOME, 'testdata', 'gr4j_cemaneige', 'pr.nc')),
             'tas': '{0}'.format(os.path.join(TESTS_HOME, 'testdata', 'gr4j_cemaneige', 'tas.nc')),
             'evap': '{0}'.format(os.path.join(TESTS_HOME, 'testdata', 'gr4j_cemaneige', 'evap.nc'))},
            'raven-gr4j-cemaneige': '{0}'.format(os.path.join(TESTS_HOME, 'testdata', 'raven-gr4j-cemaneige',
                                                              'Salmon-River-Near-Prince-George_meteo_daily.nc'))
            }
TESTDATA['raven-hmets'] = TESTDATA['raven-gr4j-cemaneige']


class WpsTestClient(WpsClient):

    def get(self, *args, **kwargs):
        query = "?"
        for key, value in kwargs.items():
            query += "{0}={1}&".format(key, value)
        return super(WpsTestClient, self).get(query)


def client_for(service):
    return WpsTestClient(service, WpsTestResponse)


def synthetic_gr4j_inputs(path):
    time = pd.date_range(start='2000-07-01', end='2002-07-01', freq='D')

    pr = 3 * np.ones(len(time))
    pr = xr.DataArray(pr, coords={'time': time}, dims='time', name='pr')
    pr.to_netcdf(os.path.join(path, 'pr.nc'))

    tas = 280 + 20 * np.cos(np.arange(len(time)) * 2 * np.pi / 365.)
    tas = xr.DataArray(tas, coords={'time': time}, dims='time', name='tas')
    tas.to_netcdf(os.path.join(path, 'tas.nc'))

    evap = 3 + 3 * np.cos(-30 + np.arange(len(time)) * 2 * np.pi / 365.)
    evap = xr.DataArray(evap, coords={'time': time}, dims='time', name='evap')
    evap.to_netcdf(os.path.join(path, 'evap.nc'))
