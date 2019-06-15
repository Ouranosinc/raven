from pywps import Service
from pywps.tests import assert_response_success

from raven.processes import TSStatsProcess, FreqAnalysisProcess, FitProcess, BaseFlowIndexProcess
from .common import client_for, TESTDATA, CFG_FILE, get_output
import xarray as xr


def test_tsstats_process():
    client = client_for(Service(processes=[TSStatsProcess(), ], cfgfiles=CFG_FILE))

    datainputs = "da=files@xlink:href=file://{da};" \
                 "freq={freq};" \
                 "op={op};" \
                 "season={season};" \
                 "variable={v};" \
        .format(da=TESTDATA['simfile_single'], freq='YS', op='max', season='JJA', v='q_sim')

    resp = client.get(
        service='WPS', request='Execute', version='1.0.0', identifier='ts_stats',
        datainputs=datainputs)

    assert_response_success(resp)
    out = get_output(resp.xml)['output']
    xr.open_dataset(out[7:])


def test_freqanalysis_process():
    client = client_for(Service(processes=[FreqAnalysisProcess(), ], cfgfiles=CFG_FILE))

    datainputs = "da=files@xlink:href=file://{da};" \
                 "t={t1};" \
                 "t={t2};" \
                 "dist={dist};" \
                 "freq={freq};" \
                 "mode={mode};" \
                 "season={season};" \
                 "variable={v};" \
        .format(da=TESTDATA['simfile_single'], freq='YS', mode='max', t1=2, t2=50, dist="gumbel_r", season='JJA',
                v='q_sim')

    resp = client.get(
        service='WPS', request='Execute', version='1.0.0', identifier='freq_analysis',
        datainputs=datainputs)

    assert_response_success(resp)
    out = get_output(resp.xml)['output']
    ds = xr.open_dataset(out[7:])
    assert ds.freq_analysis.shape == (2, 1)


def test_fit_process(ts_stats):
    client = client_for(Service(processes=[FitProcess(), ], cfgfiles=CFG_FILE))

    datainputs = "da=files@xlink:href=file://{da};" \
                 "dist={dist};" \
        .format(da=ts_stats, dist="gumbel_r")

    resp = client.get(
        service='WPS', request='Execute', version='1.0.0', identifier='fit',
        datainputs=datainputs)

    assert_response_success(resp)
    out = get_output(resp.xml)['output']
    ds = xr.open_dataset(out[7:])
    assert ds.params.shape == (2, 1)


def test_baseflowindex_process():
    client = client_for(Service(processes=[BaseFlowIndexProcess(), ], cfgfiles=CFG_FILE))

    datainputs = "q=files@xlink:href=file://{da};" \
                 "freq={freq};" \
                 "variable={v};" \
        .format(da=TESTDATA['simfile_single'], freq='YS', v='q_sim')

    resp = client.get(
        service='WPS', request='Execute', version='1.0.0', identifier='base_flow_index',
        datainputs=datainputs)

    assert_response_success(resp)
    out = get_output(resp.xml)['output']
    xr.open_dataset(out[7:])
