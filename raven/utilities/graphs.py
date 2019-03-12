from matplotlib import pyplot as plt
import xarray as xr


def ensemble_uncertainty(fns):
    """Create a graphic of the mean hydrological cycle for each model simulation.

    Parameters
    ----------
    fns : sequence
      Files storing raven hydrographs.

    Returns
    -------
    fig : plt.Figure
      Matplotlib figure showing the mean annual cycle for each simulation.
    """
    ds = [xr.open_dataset(fn) for fn in fns]
    q_sims = [h.q_sim for h in ds]
    ahs = [q.groupby('time.dayofyear').mean() for q in q_sims]

    fig, ax = plt.subplots(1, 1)
    [ah.plot(ax=ax) for ah in ahs]

    return fig
