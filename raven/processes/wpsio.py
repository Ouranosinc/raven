"""
This module contains the WPS inputs and outputs that are reused across multiple WPS processes.


"""

from pywps import LiteralInput, LiteralOutput, ComplexInput, ComplexOutput
from pywps import FORMATS, Format

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
                  max_occurs=100,
                  supported_formats=[FORMATS.NETCDF, FORMATS.TEXT, FORMATS.SHP, ])

conf = ComplexInput('conf', 'Raven/Ostrich configuration files',
                    abstract="Model configuration files, including the primary input file (rvi), the parameter "
                             "input file (rvp), the basin definition file (rvh), the time series input file "
                             "(rvt), the initial conditions file (rvc). For Ostrich, include the Ostrich "
                             "calibration config (txt) and templates (tpl).",
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

start_date = LiteralInput('start_date', 'Simulation start date (AAAA-MM-DD)',
                          abstract='Start date of the simulation (AAAA-MM-DD). '
                                   'Defaults to the start of the forcing file. ',
                          data_type='dateTime',
                          default='0001-01-01 00:00:00',
                          min_occurs=0)

end_date = LiteralInput('end_date', 'Simulation end date (AAAA-MM-DD)',
                        abstract='End date of the simulation (AAAA-MM-DD). '
                                 'Defaults to the end of the forcing file.',
                        data_type='dateTime',
                        default='0001-01-01 00:00:00',
                        min_occurs=0)

duration = LiteralInput('duration', 'Simulation duration (days)',
                        abstract='Number of simulated days, defaults to the length of the input forcings.',
                        data_type='nonNegativeInteger',
                        default=0,
                        min_occurs=0)

run_name = LiteralInput('run_name', 'Simulation name',
                        abstract='The name given to the simulation, for example <watershed>_<experiment>',
                        data_type='string',
                        default='raven-gr4j-cemaneige-sim',
                        min_occurs=0)

name = LiteralInput('name', 'Watershed name',
                    abstract='The name of the watershed the model is run for.',
                    data_type='string',
                    default='watershed',
                    min_occurs=0)

area = LiteralInput('area', 'Watershed area (km2)',
                    abstract='Watershed area (km2)',
                    data_type='float',
                    default=0.,
                    min_occurs=0)

latitude = LiteralInput('latitude', 'Latitude',
                        abstract="Watershed's centroid latitude",
                        data_type='float',
                        min_occurs=1)

longitude = LiteralInput('longitude', 'Longitude',
                         abstract="Watershed's centroid longitude",
                         data_type='float',
                         min_occurs=1)

elevation = LiteralInput('elevation', 'Elevation (m)',
                         abstract="Watershed's mean elevation (m)",
                         data_type='float',
                         min_occurs=1)

model_name = LiteralInput('model_name', 'Hydrological model identifier',
                          abstract="Hydrological model identifier: {HMETS, GR4JCN, MOHYSE}",
                          data_type='string',
                          allowed_values=('HMETS', 'GR4JCN', 'MOHYSE'),
                          min_occurs=1)


# --- #

hydrograph = ComplexOutput('hydrograph', 'Hydrograph time series (mm)',
                           supported_formats=[FORMATS.NETCDF,
                                              Format('application/zip', extension='.zip', encoding='base64')],
                           abstract='A netCDF file containing the outflow hydrographs (in m3/s) for all subbasins'
                                    'specified as `gauged` in the .rvh file. It reports period-ending time-'
                                    'averaged flows for the preceding time step, as is consistent with most '
                                    'measured stream gauge data (again, the initial flow conditions at the '
                                    'start of the first time step are included). If observed hydrographs are '
                                    'specified, they will be output adjacent to the corresponding modelled  '
                                    'hydrograph. ',
                           as_reference=True)

# TODO: adjust abstract
ensemble = ComplexOutput('ensemble', 'Multiple hydrograph time series (mm)',
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
                        supported_formats=[FORMATS.NETCDF,
                                           Format('application/zip', extension='.zip', encoding='base64')],
                        as_reference=True)

solution = ComplexOutput('solution', 'solution.rvc file to restart another simulation with the conditions '
                                     'at the end of this simulation.',
                         supported_formats=[FORMATS.TEXT,
                                            Format('application/zip', extension='.zip', encoding='base64')],
                         as_reference=True)

diagnostics = ComplexOutput('diagnostics', 'Performance diagnostic values',
                            abstract="Model diagnostic CSV file.",
                            supported_formats=[FORMATS.TEXT,
                                               Format('application/zip', extension='.zip', encoding='base64')],
                            as_reference=True)

# --- OSTRICH --- #

algorithm = LiteralInput('algorithm', 'OSTRICH Algorithm to use to calibrate model parameters',
                         abstract='Optimization algorithm to implement for this calibration run',
                         data_type='string',
                         default='DDS',
                         allowed_values=('DDS', 'SCEUA'),
                         min_occurs=0)

max_iterations = LiteralInput('max_iterations', 'Maximum number of model evaluations for the calibration run (budget)',
                              abstract='Maximum number of times OSTRICH can call the hydrological model during the '
                                       'model parameter calibrationn',
                              data_type='integer',
                              default=50,
                              allowed_values=list(range(25001)),
                              min_occurs=0)

random_seed = LiteralInput('random_seed', 'Seed for random number generator',
                           abstract="Set this value to obtain replicable results. Set to -1 to let it be random.",
                           data_type='integer',
                           default=-1,
                           min_occurs=0)

calibration = ComplexOutput('calibration', 'Ostrich calibration output',
                            abstract="Output file from Ostrich calibration run.",
                            supported_formats=[FORMATS.TEXT],
                            as_reference=True)

CalibrationResults = ComplexOutput('CalibrationResults',
                                   'ObjectiveFunction and calibrated parameters computed by Ostrich',
                                   abstract="Objective Function value after calibration using user-selected "
                                            "function, as well as the calibrated parameter set",
                                   supported_formats=[FORMATS.TEXT],
                                   as_reference=True)

calibrated_params = ComplexOutput('calibrated_params', 'Calibrated parameters',
                                  abstract="Model parameters estimated by minimizing the objective function.",
                                  supported_formats=[FORMATS.TEXT],
                                  as_reference=False)

# TODO: Add configuration files to output
# config = ComplexOutput('config', 'Configuration files',
#                        abstract="Link to configuration files.",
#                        supported_formats=)
