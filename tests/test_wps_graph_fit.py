from pywps import Service
from pywps.tests import assert_response_success

from raven.processes import GraphFitProcess
from .common import client_for, TESTDATA, CFG_FILE, get_output


def test_graph_fit(ts_stats, params):
    client = client_for(Service(processes=[GraphFitProcess(), ], cfgfiles=CFG_FILE))

    datainputs = "ts=files@xlink:href=file://{ts};" \
                 "params=files@xlink:href=file://{params};" \
        .format(ts=ts_stats, params=params)

    resp = client.get(
        service='WPS', request='Execute', version='1.0.0', identifier='graph_fit',
        datainputs=datainputs)

    assert_response_success(resp)
    out = get_output(resp.xml)
    assert out['graph_fit'].endswith('.png')
