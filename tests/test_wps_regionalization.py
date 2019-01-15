# -*- coding: utf-8 -*-

import datetime as dt
from pywps import Service
from pywps.tests import assert_response_success
from raven.processes import RegionalisationProcess
from .common import client_for, TESTDATA, get_output, urlretrieve, CFG_FILE
import xarray as xr

datainputs = "ts=files@xlink:href=file://{ts};" \
             "start_date={start_date};" \
             "end_date={end_date};" \
             "name={name};" \
             "area={area};" \
             "latitude={latitude};" \
             "longitude={longitude};" \
             "elevation={elevation};" \
             "model_name={model_name};" \
             "min_NSE={min_NSE};" \
             "ndonors={ndonors};" \
             "method={method};"

inputs = dict(start_date=dt.datetime(2000, 1, 1),
              end_date=dt.datetime(2002, 1, 1),
              name='Salmon',
              run_name='test',
              area='4250.6',
              elevation='843.0',
              latitude=40.4848,
              longitude=-103.3659,
              min_NSE=0.6,
              )


class TestRegionalisation:

    def testRegionalisationHMETS(self):

        client = client_for(Service(processes=[RegionalisationProcess(), ], cfgfiles=CFG_FILE))

        inp = inputs.copy()
        inp['ts'] = TESTDATA['raven-hmets-nc-ts']
        inp['model_name'] = 'HMETS'
        inp['ndonors'] = 4
        inp['method'] = 'MLR'

        resp = client.get(service='WPS', request='execute', version='1.0.0',
                          identifier='RegionalisationProcess',
                          datainputs=datainputs.format(**inp))

        assert_response_success(resp)
        out = get_output(resp.xml)
        assert 'hydrograph' in out
        hydrograph, _ = urlretrieve(out['hydrograph'])
        ensemble, _ = urlretrieve(out['ensemble'])

        with xr.open_dataset(hydrograph) as ds:
            assert 'Qsim' in ds

        with xr.open_dataset(ensemble) as ds:
            assert 'ens' in ds.Qsim.dims
