import xarray as xr
from raven.models import get_model
import numpy as np


def realization(n):
    """Return a realization coordinate.

    Parameters
    ----------
    n : int
      Size of the ensemble.
    """
    return xr.IndexVariable('realization',
                            data=range(n),
                            attrs={'standard_name': 'realization', 'axis': 'E', 'units': 1,
                                   'long_name': 'Label identifying the ensemble member'})


def param(model):
    """Return a parameter coordinate.

    Parameters
    ----------
    model : str
      Model name.
    """
    model = get_model(model)
    return xr.IndexVariable('param',
                            data=np.array(model.RVP.params._fields),
                            attrs={'standard_name': 'parameter', 'long_name': '{} model parameter name'.format(model)})
