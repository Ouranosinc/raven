import xarray as xr
import xclim.sdba as sdba
from xclim import subset
from . common import TESTDATA
import datetime as dt
from raven.utilities import forecasting


class TestHindcasting:
    def test_hindcasting_GEPS(self):

        shape=TESTDATA['watershed_vector']
        
        date=dt.datetime(2017,6,2)
        data=forecasting.get_hindcast_day(shape,date)