from raven.utilities import regionalization as reg
import datetime as dt
from .common import TESTDATA


def test_regionalization():
    model = 'GR4JCN'
    nash, params = reg.read_gauged_params(model)
    variables = ['latitude', 'longitude', 'area','forest']
    props = reg.read_gauged_properties(variables)
    ungauged_props = {'latitude': 40.4848, 'longitude': -103.3659,'area': 4250.6, 'forest': 0.4}

    qsim, ens = reg.regionalize('SP_IDW', model, nash, params,
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

    assert(qsim.max() > 1)
    assert(len(ens) == 2)
    assert 'realization' in ens.dims
    assert 'param' in ens.dims
