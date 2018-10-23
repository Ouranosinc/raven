from raven import models
import pandas as pd
import numpy as np
import xarray as xr
from .common import TESTDATA

class TestGR4JCemaneige():

    def test_simple(self):
        time = pd.date_range('2000-01-01', freq='D', periods=365*3)
        n = len(time)
        temp = pd.Series(-20*np.sin(np.arange(n)/365.25 * 2 * np.pi)+5, index=time)
        precip = pd.Series(((np.random.rand(n)-.5)*20).clip(0), index=time)
        evap = (temp/400 + (np.random.rand(n)-.5)*8).clip(0)
        df = pd.DataFrame({'Temp': temp, 'Evap': evap, 'Prec': precip})

        params = [300, -1, 200, 2, .5, 8]
        out = models.gr4j(df, params)
        return out


    def test_nc_data(self):
        paths = TESTDATA['gr4j-cemaneige'].values()
        ds = xr.open_mfdataset(paths).rename({'pr': 'Prec', 'tas': 'Temp', 'evap': 'Evap'})
        ds['Temp'] -= 273.15
        df = ds.to_dataframe()
        params = [300, -1, 200, 2, .5, 8]
        q = models.gr4j(df, params)
        return q
