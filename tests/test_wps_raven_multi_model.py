import datetime as dt

from pywps import Service
from pywps.tests import assert_response_success
from raven.processes import RavenMultiModelProcess, GraphEnsUncertaintyProcess, GraphSingleHydrographProcess
from .common import client_for, TESTDATA, CFG_FILE, get_output

import pytest


# @pytest.mark.dependency()
def test_raven_multi_model_process(request):
    client = client_for(Service(processes=[RavenMultiModelProcess(), ], cfgfiles=CFG_FILE))

    gr4jcn = '0.529, -3.396, 407.29, 1.072, 16.9, 0.947'
    hmets = '9.5019, 0.2774, 6.3942, 0.6884, 1.2875, 5.4134, 2.3641, 0.0973, 0.0464, 0.1998, 0.0222, -1.0919, ' \
            '2.6851, 0.3740, 1.0000, 0.4739, 0.0114, 0.0243, 0.0069, 310.7211, 916.1947'

    datainputs = "ts=files@xlink:href=file://{ts};" \
                 "gr4jcn={gr4jcn};" \
                 "hmets={hmets};" \
                 "start_date={start_date};" \
                 "end_date={end_date};" \
                 "name={name};" \
                 "run_name={run_name};" \
                 "area={area};" \
                 "latitude={latitude};" \
                 "longitude={longitude};" \
                 "elevation={elevation};" \
        .format(ts=TESTDATA['raven-gr4j-cemaneige-nc-ts'],
                gr4jcn=gr4jcn,
                hmets=hmets,
                start_date=dt.datetime(2000, 1, 1),
                end_date=dt.datetime(2002, 1, 1),
                name='Salmon',
                run_name='test',
                area='4250.6',
                elevation='843.0',
                latitude=54.4848,
                longitude=-123.3659,
                )

    resp = client.get(
        service='WPS', request='Execute', version='1.0.0', identifier='raven-multi-model',
        datainputs=datainputs)

    assert_response_success(resp)
    out = get_output(resp.xml)
    assert out['hydrograph'].endswith('.zip')
    request.config.cache.set('zipfn', out['hydrograph'])


# @pytest.mark.dependency(depends=['test_raven_multi_model_process'])
@pytest.mark.skip
def test_graph_ensemble_uncertainty(request):
    client = client_for(Service(processes=[GraphEnsUncertaintyProcess(), ], cfgfiles=CFG_FILE))

    zipfn = request.config.cache.get('zipfn', None)

    resp = client.get(
        service='WPS', request='Execute', version='1.0.0', identifier='graph_ensemble_uncertainty',
        datainputs="sims=files@xlink:href=file://{};".format(zipfn))

    assert_response_success(resp)


@pytest.mark.skip
def test_graph_single_hydrograph(request):
    client = client_for(Service(processes=[GraphSingleHydrographProcess(), ], cfgfiles=CFG_FILE))
    datainputs = "sim=files@xlink:href=file://{sim};".format(sim=TESTDATA['simfile_single'])

    resp = client.get(
        service='WPS', request='Execute', version='1.0.0', identifier='graph_single_hydrograph',
        datainputs=datainputs)

    assert_response_success(resp)
