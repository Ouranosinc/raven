import numpy as np
import xarray as xr
from ravenpy.utilities import graphs
from ravenpy.utilities.testdata import get_local_testdata
from xclim.indicators.land import fit, stats


def test_ts_fit_graph():
    fn = get_local_testdata(
        "hydro_simulations/raven-gr4j-cemaneige-sim_hmets-0_Hydrographs.nc"
    )
    ds = xr.open_dataset(fn)

    ts = stats(ds.q_sim, op="max", freq="M")
    p = fit(ts)
    np.testing.assert_array_equal(p.isnull(), False)

    fig = graphs.ts_fit_graph(ts, p)
    return fig
