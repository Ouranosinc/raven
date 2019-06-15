from importlib import reload
from raven.utilities import graphs
reload(graphs)
from .common import TESTDATA
from raven.processes.wps_q_stats import fit, stats
import xarray as xr


def test_ts_fit_graph():
    fn = TESTDATA['simfile_single']
    ds = xr.open_dataset(fn)

    ts = stats(ds.q_sim, op='max')
    p = fit(ts)

    fig = graphs.ts_fit_graph(ts, p)
    return fig
