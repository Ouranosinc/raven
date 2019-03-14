from matplotlib import pyplot as plt
import xarray as xr

#TODO: Ensure that the nameList items are in the same order as the generator netCDF files!!!
    
def ensemble_uncertainty_annual(fns,nameList):
    """Create a graphic of the mean hydrological cycle for each model simulation.

    Parameters
    ----------
    fns : sequence
      Files storing raven hydrographs.
    nameList: list
      List of strings of the model runs to add to the legend.
    Returns
    -------
    fig : plt.Figure
      Matplotlib figure showing the mean annual cycle for each simulation.
    """
    ds = [xr.open_dataset(fn) for fn in fns]
    simulated_streamflow = [h.q_sim for h in ds]
    ahs = [q.groupby('time.dayofyear').mean() for q in simulated_streamflow]

    fig, ax = plt.subplots(1, 1)
    [ah.plot(ax=ax) for ah in ahs]
    ax.legend(nameList)

    return fig
