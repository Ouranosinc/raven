from ravenpy.utilities import graphs
from .common import TESTDATA
from xclim.indicators.land import fit, stats
import xarray as xr


def test_ts_fit_graph():
    fn = TESTDATA['simfile_single']
    ds = xr.open_dataset(fn)

    ts = stats(ds.q_sim, op='max').dropna("time")
    p = fit(ts)

    fig = graphs.ts_fit_graph(ts, p)
    return fig
