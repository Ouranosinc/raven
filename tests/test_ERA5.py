import pytest
import datetime as dt
import numpy as np
import xarray as xr

from pywps import Service
from pywps.tests import assert_response_success

from . common import client_for, TESTDATA, CFG_FILE, get_output, urlretrieve
from raven.processes import RavenHMETSProcess
from raven.models import HMETS
import json
import matplotlib.pyplot as plt

params = (9.5019, 0.2774, 6.3942, 0.6884, 1.2875, 5.4134, 2.3641, 0.0973, 0.0464, 0.1998, 0.0222, -1.0919,
          2.6851, 0.3740, 1.0000, 0.4739, 0.0114, 0.0243, 0.0069, 310.7211, 916.1947)


@pytest.mark.slow
class TestRavenERA5:
    def test_simple(self, era5_hr):
        model = HMETS()
        model(ts=era5_hr,
              params=params,
              start_date=dt.datetime(2018, 1, 1),
              end_date=dt.datetime(2018, 8, 10),
              name='Salmon',
              run_name='test-hmets-era5',
              area=4250.6,
              elevation=843.0,
              latitude=54.4848,
              longitude=-123.3659,
              rain_snow_fraction="RAINSNOW_DINGMAN",
              tas={'linear_transform': (1.0, -273.15), 'time_shift': -.25},
              pr={'linear_transform': (24000.0, 0.0), 'time_shift': -.25}
              )


@pytest.mark.slow
class TestRavenERA5Process:

    def test_simple(self, era5_hr):
        client = client_for(Service(processes=[RavenHMETSProcess(), ], cfgfiles=CFG_FILE))

        datainputs = "ts=files@xlink:href=file://{ts};" \
                     "params={params};" \
                     "start_date={start_date};" \
                     "end_date={end_date};" \
                     "init={init};" \
                     "name={name};" \
                     "run_name={run_name};" \
                     "area={area};" \
                     "latitude={latitude};" \
                     "longitude={longitude};" \
                     "elevation={elevation};" \
                     "rain_snow_fraction={rain_snow_fraction};"\
                     "nc_spec={tas};"\
                     "nc_spec={pr}"\
            .format(ts=era5_hr,
                    params=params,
                    start_date=dt.datetime(2018, 1, 1),
                    end_date=dt.datetime(2018, 8, 10),
                    init='155,455',
                    name='Salmon',
                    run_name='test-hmets-era5',
                    area='4250.6',
                    elevation='843.0',
                    latitude=54.4848,
                    longitude=-123.3659,
                    rain_snow_fraction="RAINSNOW_DINGMAN",
                    pr=json.dumps({'pr': {'linear_transform': (24000.0, 0.0), 'time_shift': -.25}}),
                    tas=json.dumps({'tas': {'linear_transform': (1.0, -273.15), 'time_shift': -.25}}),
                    )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='raven-hmets',
            datainputs=datainputs)
        assert_response_success(resp)
#        out = get_output(resp.xml)
#        tmp,_= urlretrieve(out['hydrograph'])
#        q=xr.open_dataset(tmp)
#        plt.plot(q['q_sim'])
#        plt.show()
