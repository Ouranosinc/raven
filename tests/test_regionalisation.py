from raven.utilities import regionalization as reg
import datetime as dt
from .common import TESTDATA


def test_regionalization():
    model = 'GR4JCN'
    nash, params = reg.read_gauged_params(model)
    variables = ['longitude', 'latitude']
    props = reg.read_gauged_properties()[variables]
    ungauged_props = {'longitude': .7, 'latitude': .7}

    qsim, ensemble = reg.regionalize('SP_IDW', model, nash, params,
                    props, ungauged_props,
                    start_date=dt.datetime(2000, 1, 1),
                    end_date=dt.datetime(2002, 1, 1),
                    name='Salmon',
                    run_name='test',
                    area='4250.6',
                    elevation='843.0',
                    latitude=40.4848,
                    longitude=-103.3659,
                    min_NSE=0.6,
                    size=2,
                    ts=TESTDATA['raven-hmets-nc-ts'])
    
    assert(qsim.q_sim.max()>1)
    assert(len(ensemble.q_sim)==2)
