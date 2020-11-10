from . common import TESTDATA
import datetime as dt
from raven.utilities import forecasting
from raven.models import GR4JCN
import pdb

class TestHindcasting:
    def test_hindcasting_GEPS(self):

        shape=TESTDATA['watershed_vector']
        
        date=dt.datetime(2017,6,1)
        fcst=forecasting.get_hindcast_day(shape,date,climate_model='GEPS')
        fcst.to_netcdf('/home/ets/src/raventest/raven/tests/hindcastfile.nc')
        
        ts = TESTDATA["raven-gr4j-cemaneige-nc-ts"]
        model = GR4JCN()
        model(
            ts,
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2005, 6, 1),
            area=44250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            params=(0.529, -3.396, 407.29, 1.072, 16.9, 0.947),
        )
        rvc=model.outputs["solution"]
        model.resume(rvc)
        
        model(ts=('/home/ets/src/raventest/raven/tests/hindcastfile.nc'),start_date=dt.datetime(2017, 6, 1),
            end_date=dt.datetime(2017, 6, 2),
            area=44250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            params=(0.529, -3.396, 407.29, 1.072, 16.9, 0.947),overwrite=True)
        pdb.set_trace()