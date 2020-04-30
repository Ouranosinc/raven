import pytest
import datetime as dt
import numpy as np
import xarray as xr
import json
import matplotlib.pyplot as plt
from pywps import Service
from pywps.tests import assert_response_success

import tempfile
from . common import client_for, TESTDATA, CFG_FILE, get_output, urlretrieve
from raven.processes import RavenHMETSProcess
from raven.models import HMETS


NRCAN_path='https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/birdhouse/1-Datasets/gridded_obs/nrcan_v2.ncml'

# Temporary path
filepath = tempfile.mkdtemp() + "/NRCAN_ts.nc"

# Get information for given catchment, could be passed in parameter to the function
salmon=xr.open_dataset(TESTDATA['raven-hmets-nc-ts'])
lat = salmon.lat.values[0]
lon = salmon.lon.values[0]

#Start and end dates
start_date=dt.datetime(2006, 1, 1)
end_date=dt.datetime(2007, 12, 31)
              
# Get data for the covered period including a 1-degree bounding box for good measure. Eventually we will be
# able to take the catchment polygon as a mask and average points residing inside.

ds=xr.open_dataset(NRCAN_path).sel(lat=slice(lat+1,lat-1), lon=slice(lon-1,lon+1), time=slice(start_date, end_date)).mean(dim={'lat','lon'}, keep_attrs=True)
ds.to_netcdf(filepath)


class TestRavenNRCAN:
    def test_simple(self):

        
        params = (9.5019, 0.2774, 6.3942, 0.6884, 1.2875, 5.4134, 2.3641, 0.0973, 0.0464, 0.1998, 0.0222, -1.0919,
          2.6851, 0.3740, 1.0000, 0.4739, 0.0114, 0.0243, 0.0069, 310.7211, 916.1947)
           
        model = HMETS()
        model(ts=filepath,
              params=params,
              start_date=start_date,
              end_date=end_date,
              name='Salmon',
              run_name='test-hmets-NRCAN',
              area=4250.6,
              elevation=843.0,
              latitude=54.4848,
              longitude=-123.3659,
              rain_snow_fraction="RAINSNOW_DINGMAN",
              tasmax={'linear_transform': (1.0, -273.15)},
              tasmin={'linear_transform': (1.0, -273.15)},
              pr={'linear_transform': (86400, 0.0)}
              )
        
        
class TestRavenNRCANProcess:

    def test_simple(self):
        client = client_for(Service(processes=[RavenHMETSProcess(), ], cfgfiles=CFG_FILE))

        params = '9.5019, 0.2774, 6.3942, 0.6884, 1.2875, 5.4134, 2.3641, 0.0973, 0.0464, 0.1998, 0.0222, -1.0919, ' \
                 '2.6851, 0.3740, 1.0000, 0.4739, 0.0114, 0.0243, 0.0069, 310.7211, 916.1947'
     
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
                     "nc_spec={tasmax};"\
                     "nc_spec={tasmin};"\
                     "nc_spec={pr}"\
            .format(ts=filepath, # This is a file on disk, need to pass 'g'
                    params=params,
                    start_date=start_date,
                    end_date=end_date,
                    init='155,455',
                    name='Salmon',
                    run_name='test-hmets-NRCAN',
                    area='4250.6',
                    elevation='843.0',
                    latitude=lat,
                    longitude=lon,
                    rain_snow_fraction="RAINSNOW_DINGMAN",
                    tasmax=json.dumps({'tasmax': {'linear_transform': (1.0, -273.15)}}),
                    tasmin=json.dumps({'tasmin': {'linear_transform': (1.0, -273.15)}}),
                    pr=json.dumps({'pr': {'linear_transform': (86400.0, 0.0)}}),
                    )
            
        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='raven-hmets',
            datainputs=datainputs)

        assert_response_success(resp)
        
        # For testing purposes
        #out = get_output(resp.xml)
        #tmp_file, _ = urlretrieve(out['hydrograph'])
        #tmp=xr.open_dataset(tmp_file)
        #plt.plot(tmp['q_sim'])
        #plt.show()
        
