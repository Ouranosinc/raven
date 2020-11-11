from . common import TESTDATA
import datetime as dt
from raven.utilities import forecasting
from raven.models import GR4JCN


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
            end_date=dt.datetime(2002, 6, 1),
            area=44250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            params=(0.529, -3.396, 407.29, 1.072, 16.9, 0.947),
        )
        rvc=model.outputs["solution"]       
        
        # It is necessary to clean the model state because the input variables of the previous
        # model are not the same as the ones provided in the forecast model. therefore, if we
        # do not clean, the model will simply add the hindcast file to the list of available
        # data provided in the testdata above. Then the dates will not work, and the model errors.
        model = GR4JCN()
        
        # Now, relaunch the model at the good states.
        model.rvc.parse(rvc.read_text())
        
        model(ts=('/home/ets/src/raventest/raven/tests/hindcastfile.nc'),start_date=dt.datetime(2017, 6, 1),
            end_date=dt.datetime(2017, 6, 2),
            area=44250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            params=(0.529, -3.396, 407.29, 1.072, 16.9, 0.947),overwrite=True)
        
        assert len(model.q_sim.values)==2
        