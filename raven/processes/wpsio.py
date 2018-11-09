"""
This module contains the WPS inputs and outputs that are reused across multiple WPS processes.


"""

from pywps import LiteralInput, LiteralOutput, ComplexInput, ComplexOutput
from pywps import FORMATS

# ---------------------------------------- #
# ---------------- Inputs ---------------- #
# ---------------------------------------- #


ts = ComplexInput('ts', 'Input time series files',
                  abstract='Files (text or netCDF) storing'
                           'daily liquid precipitation (pr), '
                           'solid precipitation (prsn), '
                           'minimum temperature (tasmin), '
                           'maximum temperature (tasmax), '
                           'potential evapotranspiration (evspsbl) and '
                           'observed streamflow (qobs [m3/s]).',
                  min_occurs=1,
                  supported_formats=[FORMATS.NETCDF, FORMATS.TEXT])

conf = ComplexInput('conf', 'Raven configuration files',
                    abstract="Model configuration files, including the primary input file (rvi), the parameter "
                             "input file (rvp), the basin definition file (rvh), the time series input file "
                             "(rvt), and the initial conditions file (rvc). ",
                    min_occurs=5,
                    max_occurs=5,
                    supported_formats=[FORMATS.TEXT])

rvi = ComplexInput('rvi', 'Primary input file',
                   abstract="The primary input file stores the model simulation options and numerical options.",
                   min_occurs=1,
                   max_occurs=1,
                   supported_formats=[FORMATS.TEXT])

rvp = ComplexInput('rvp', 'Classed parameter input file',
                   abstract="The classed parameter input file stores a database of soil, vegetation, river, "
                            "aquifer, and land class pro-perties. Not all classes specified in the *.rvp file "
                            "need to be included in the model.",
                   min_occurs=1,
                   max_occurs=1,
                   supported_formats=[FORMATS.TEXT])

rvh = ComplexInput('rvh', 'HRU / Basin definition file',
                   abstract="The HRU/basin definition file describes the topology of the basin network and the "
                            "class membership of all constituent HRUs.",
                   min_occurs=1,
                   max_occurs=1,
                   supported_formats=[FORMATS.TEXT])

rvt = ComplexInput('rvt', 'Time series input file',
                   abstract="The time series input file is used to store time series of forcing functions ("
                            "precipitation, temperature, etc.).",
                   min_occurs=1,
                   max_occurs=1,
                   supported_formats=[FORMATS.TEXT])

rvc = ComplexInput('rvc', 'Initial conditions input file',
                   abstract="The initial conditions input file is used to store the initial conditions for the "
                            "model. By default, the initial conditions for all model state variables is zero, "
                            "and there are no required commands in this file (it could even be completely "
                            "empty).",
                   min_occurs=0,
                   default="",
                   supported_formats=[FORMATS.TEXT])

# ---------------------------------------- #
# ---------------- Outputs --------------- #
# ---------------------------------------- #

hydrograph = ComplexOutput('hydrograph', 'Hydrograph time series (mm)',
                           supported_formats=[FORMATS.NETCDF],
                           abstract='A netCDF file containing the outflow hydrographs (in m3/s) for all subbasins'
                                    'specified as `gauged` in the .rvh file. It reports period-ending time-'
                                    'averaged flows for the preceding time step, as is consistent with most '
                                    'measured stream gauge data (again, the initial flow conditions at the '
                                    'start of the first time step are included). If observed hydrographs are '
                                    'specified, they will be output adjacent to the corresponding modelled  '
                                    'hydrograph. ',
                           as_reference=True)

storage = ComplexOutput('storage', 'Watershed storage time series (mm)',
                        abstract='A netCDF file describing the total storage of water (in mm) in all water '
                                 'storage compartments for each time step of the simulation. Mass balance '
                                 'errors, cumulative input (precipitation), and output (channel losses) are '
                                 'also included. Note that the precipitation rates in this file are '
                                 'period-ending, i.e., this is the precipitation rate for the time step '
                                 'preceding the time stamp; all water storage variables represent '
                                 'instantaneous reports of the storage at the time stamp indicate.',
                        supported_formats=[FORMATS.NETCDF],
                        as_reference=True)

solution = ComplexOutput('solution', 'solution.rvc file to restart another simulation with the conditions '
                                     'at the end of this simulation.',
                         supported_formats=[FORMATS.TEXT],
                         as_reference=True)

diagnostics = ComplexOutput('diagnostics', 'Performance diagnostic values',
                            abstract="Model diagnostic CSV file.",
                            supported_formats=[FORMATS.TEXT],
                            as_reference=True)
