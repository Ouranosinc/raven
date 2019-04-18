from pywps import Service
from pywps.tests import assert_response_success

from raven.processes import TSStatsProcess
from .common import client_for, TESTDATA, CFG_FILE, get_output
import xarray as xr


def test_tsstats_process():

    client = client_for(Service(processes=[TSStatsProcess(), ], cfgfiles=CFG_FILE))

    datainputs = "da=files@xlink:href=file://{da};"\
                 "freq={freq};"\
                 "op={op};"\
                 "season={season};"\
                 "variable={v};"\
                 .format(da=TESTDATA['simfile_single'], freq='YS', op='max', season='JJA', v='q_sim')

    resp = client.get(
        service='WPS', request='Execute', version='1.0.0', identifier='ts_stats',
        datainputs=datainputs)

    assert_response_success(resp)
    out = get_output(resp.xml)['output']
    xr.open_dataset(out[7:])
