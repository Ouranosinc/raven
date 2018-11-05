import os
import raven
from raven.models import raven_templates
import xarray as xr
import csv

# Model executable
raven_exec = os.path.abspath(os.path.join(os.path.dirname(raven.__file__), '..', 'bin', 'raven'))

# Dictionary of potential variable names, keyed by CF standard name.
# http://cfconventions.org/Data/cf-standard-names/60/build/cf-standard-name-table.html
# PET is the potential evapotranspiration, while evspsbl is the actual evap.
variable_names = {'tasmin': ['tasmin', 'tmin'],
                  'tasmax': ['tasmax', 'tmax'],
                  'pr': ['pr', 'precip', 'prec', 'rain', 'rainfall', 'precipitation'],
                  'prsn': ['prsn', 'snow', 'snowfall', 'solid_precip'],
                  'evspsbl': ['pet', 'evap', 'evapotranspiration'],
                  'water_volume_transport_in_river_channel': ['qobs', 'discharge', 'streamflow']
                  }


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
      The command line arguments to launch the model with the configuration files.
    """
    inpath = raven_templates[name]
    model_path = os.path.join(outpath, 'model')
    output_path = os.path.join(outpath, 'output')

    # Create subdirectory
    os.mkdir(model_path)
    os.mkdir(output_path)

    for ext, param in params.items():
        fn = name + os.path.extsep + ext
        txt = rv_format(os.path.join(inpath, fn), param)

        with open(os.path.join(model_path, fn), 'w') as f:
            f.write(txt)

    cmd = os.path.join(model_path, 'raven')
    os.symlink(raven_exec, cmd)

    return [cmd, os.path.join(outpath, 'model', name), '-o', os.path.join(outpath, 'output')]


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
    """Find for each variable the file storing it's data and the name of the netCDF variable.

    Parameters
    ----------
    fns : sequence
      Paths to netCDF files.
    variables : sequence
      Names of the variables to look for. Specify their CF standard name, a dictionary of
      alternative names will be used for the lookup.

    Returns
    -------
    out : list
      A dictionary keyed by variable storing the file storing each variable.
    """
    out = {}
    for fn in fns:
        ds = xr.open_dataset(fn)
        for var in variables:
            for alt_name in variable_names[var]:
                if alt_name in ds.data_vars:
                    out[var] = fn
                    out[var + '_var'] = alt_name
                    break

    for var in variables:
        if var not in out.keys():
            raise ValueError("{} not found in files.".format(var))

    return out


def read_diagnostics(f):
    """Return a dictionary representation of the output diagnostic file.


    Parameters
    ----------
    f : file
      Diagnostic file.
    """

    reader = csv.reader(f.readlines())
    header = reader.next()
    content = reader.next()

    out = dict(zip(header, content))
    out.pop('')
    return out
