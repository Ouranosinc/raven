import os
import raven
from raven.models import raven_templates
import xarray as xr

# Model executable
exec = os.path.join(os.path.abspath(os.path.dirname(raven.__file__)), 'bin', 'raven')

def rv_format(fn, kwds):
    """Read the model input template file and fill the given arguments."""

    if None in kwds.values():
        raise ValueError("Some parameters are not properly set.")

    with open(fn) as f:
        txt = f.read()

    return txt.format(**kwds)


def setup_model(name, outpath, params):
    """Create a subdirectory and write the RAVEN configuration files.

    Parameters
    ----------
    name : str {'raven-gr4j'}
      Model name.
    outpath : str
      Work directory where the configuration files will be written.
    params : dict
      Model parameter values to be written into template.

    Returns
    -------
    cmd : str
      The model executable symbolic link in created directory.
    """
    inpath = raven_templates[name]

    # Create subdirectory
    os.mkdir(os.path.join(outpath, 'model'))
    os.mkdir(os.path.join(outpath, 'output'))

    for ext, param in params.items():
        fn = name + os.path.extsep + ext
        txt = rv_format(os.path.join(inpath, fn), param)

        with open(os.path.join(outpath, 'model', fn), 'w') as f:
            f.write(txt)

    cmd = os.path.join(outpath, 'model', 'raven')
    os.symlink(exec, cmd)

    return cmd


def start_date(fns):
    """Return the common starting date and time of netCDF files.

    Parameters
    ----------
    fns : sequence
      Sequence of netCDF file names for forcing data.

    Returns
    -------
    start : datetime
      The start datetime of the forcing files.
    """
    ds = xr.open_mfdataset(fns)
    return ds.indexes['time'][0]
