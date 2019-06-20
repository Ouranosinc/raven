# -*- coding: utf-8 -*-

import pytest
import datetime as dt
from pywps import Service
from pywps.tests import assert_response_success
from raven.processes import RegionalisationProcess
from .common import client_for, TESTDATA, get_output, urlretrieve, CFG_FILE
import xarray as xr
import json


datainputs = "ts=files@xlink:href=file://{ts};" \
             "start_date={start_date};" \
             "end_date={end_date};" \
             "name={name};" \
             "latitude={latitude};" \
             "longitude={longitude};" \
             "model_name={model_name};" \
             "min_NSE={min_NSE};" \
             "ndonors={ndonors};" \
             "method={method};" \
             "properties={properties};" \
             "area={area};" \
             "elevation={elevation};"

inputs = dict(start_date=dt.datetime(2000, 1, 1),
              end_date=dt.datetime(2002, 1, 1),
              name='Salmon',
              run_name='test',
              latitude=45,
              longitude=-80,
              min_NSE=0.6,
              properties=json.dumps({'latitude': 45, 'longitude': -80, 'forest': 0.7}),
              area=5600,
              elevation=100,
              )


class TestRegionalisation:

    @pytest.mark.parametrize("method", ('MLR', 'SP', 'PS', 'SP_IDW', 'PS_IDW', 'SP_IDW_RA', 'PS_IDW_RA'))
    def testRegionalisationHMETS(self, method):

        client = client_for(Service(processes=[RegionalisationProcess(), ], cfgfiles=CFG_FILE))

        inp = inputs.copy()
        inp['ts'] = TESTDATA['raven-hmets-nc-ts']
        inp['model_name'] = 'GR4JCN'
        inp['ndonors'] = 2
        inp['method'] = method

        resp = client.get(service='WPS', request='execute', version='1.0.0',
                          identifier='regionalisation',
                          datainputs=datainputs.format(**inp))

        assert_response_success(resp)
        out = get_output(resp.xml)
        assert 'hydrograph' in out
        hydrograph, _ = urlretrieve(out['hydrograph'])
        ensemble, _ = urlretrieve(out['ensemble'])

        with xr.open_dataset(hydrograph) as ds:
            assert 'q_sim' in ds

        with xr.open_dataset(ensemble) as ds:
            assert 'realization' in ds.q_sim.dims

    def test_notebook(self):
        client = client_for(Service(processes=[RegionalisationProcess(), ], cfgfiles=CFG_FILE))

        inp = inputs.copy()
        inp.update(
            ts=TESTDATA['raven-hmets-nc-ts'],
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2002, 1, 1),
            area=4250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            method='PS',  # One of the methods described above
            model_name='HMETS',  # One of the two models are allowed: HMETS and GR4JCN
            min_nse=0.7,
            ndonors=5,  # Number of donors we want to use. Usually between 4 and 8 is a robust number.
            properties=json.dumps({'latitude': 54.4848, 'longitude': -123.3659, 'forest': 0.4}),
        )

        resp = client.get(service='WPS', request='execute', version='1.0.0',
                          identifier='regionalisation',
                          datainputs=datainputs.format(**inp))

        assert_response_success(resp)
