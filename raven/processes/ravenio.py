import os
import raven
from raven.models import raven_templates
import xarray as xr

# Model executable
exec = os.path.abspath(os.path.join(os.path.dirname(raven.__file__), '..', 'bin', 'raven'))

def rv_format(fn, kwds):
    """Read the model input template file and fill the given arguments."""
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


def start_end_date(fns):
    """Return the common starting and ending date and time of netCDF files.

    Parameters
    ----------
    fns : sequence
      Sequence of netCDF file names for forcing data.

    Returns
    -------
    start : datetime
      The first datetime of the forcing files.
    end : datetime
      The last datetime of the forcing files.
    """
    ds = xr.open_mfdataset(fns)
    return ds.indexes['time'][0], ds.indexes['time'][-1]


def assign_files(fns, variables):
    """Find for each variable the file storing it's data.

    Parameters
    ----------
    fns : sequence
      Paths to netCDF files.
    variables : sequence
      Names of the variables to look for.

    Returns
    -------
    out : list
      A dictionary keyed by variable storing the file storing each variable.
    """
    out = {}
    for fn in fns:
        ds = xr.open_dataset(fn)
        for var in variables:
            if var in ds:
                out[var] = fn

    for var in variables:
        if var not in out.keys():
            raise ValueError("{} not found in files.".format(var))

    return [out[v] for v in variables]



