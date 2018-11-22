import datetime as dt
import logging
import os
# from pywps.app.Common import Metadata
import subprocess
from collections import OrderedDict as Odict

from pywps import ComplexInput, ComplexOutput
from pywps import FORMATS
from pywps import LiteralInput
from pywps import Process

from . import ravenio

LOGGER = logging.getLogger("PYWPS")

"""
Notes
-----

The configuration files for RAVEN's GR4J-Cemaneige model and in models/raven-gr4j-cemaneige.
All parameters that could potentially be user-defined are tagged using {}. These tags need to be replaced by
actual values before the model is launched. There are two mechanisms to make those replacements:

1. Set-up frozen values in the `defaults` dictionary (see below). These cannot be modified by the WPS request.
2. Provide values through the WPS request (request- or default-defined).

Multiple processes running the GR4J-Cemaneige model can thus template the model differently.
"""

defaults = Odict(
    rvi=dict(run_name=None, Start_Date=None, End_Date=None, Duration=None, TimeStep=1.0,
             EvaluationMetrics='NASH_SUTCLIFFE RMSE'),
    rvp=dict(GR4J_X1=None, GR4J_X2=None, GR4J_X3=None, GR4J_X4=None, AvgAnnualSnow=None, AirSnowCoeff=None),
    rvc=dict(SOIL_0=None, SOIL_1=None),
    rvh=dict(NAME=None, AREA=None, ELEVATION=None, LATITUDE=None, LONGITUDE=None),
    rvt=dict(pr=None, prsn=None, tasmin=None, tasmax=None, evspsbl=None, water_volume_transport_in_river_channel=None)
)


class RavenGR4JCemaNeigeProcess(Process):
    """
    RAVEN emulator for the GR4J-Cemaneige model.

    This process runs the GR4J-Cemaneige model using a RAVEN emulator. Users need to provide netCDF input files for
    rain, snow minimum and maximum temperature as well as potential evapotranspiration. To run diagnostics, observed
    stream flows are also required.


    """
    identifier = 'raven-gr4j-cemaneige'
    abstract = 'Raven GR4J + CEMANEIGE hydrological model'
    title = ''
    version = ''
    defaults = defaults  # Model defaults
    pdefaults = Odict([('SOIL_PROD', 0.696),
                       ('GR4J_X2', 0.7),
                       ('GR4J_X3', 19.7),
                       ('GR4J_X4', 2.09),
                       ('AvgAnnualSnow', 123.3),
                       ('AirSnowCoeff', 0.75)])

    nc = ComplexInput('nc', 'netCDF input files',
                      abstract='NetCDF file or files storing'
                               ' daily liquid precipitation (pr), '
                               'solid precipitation (prsn), '
                               'minimum temperature (tasmin), '
                               'maximum temperature (tasmax), '
                               'potential evapotranspiration (evspsbl) and '
                               'observed streamflow (qobs [m3/s]).',
                      min_occurs=1,
                      supported_formats=[FORMATS.NETCDF])

    params = LiteralInput('params', 'Comma separated list of model parameters',
                          abstract='Parameters: ' + ', '.join(pdefaults.keys()),
                          data_type='string',
                          default=', '.join(str(p) for p in pdefaults.values()),
                          min_occurs=0)

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

    init = LiteralInput('init', 'Initial soil conditions',
                        abstract='Underground reservoir levels: SOIL_0, SOIL_1',
                        data_type='string',
                        default='0, 0',
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

    # --- #

    hydrograph = ComplexOutput('hydrograph', 'Hydrograph time series (mm)',
                               supported_formats=FORMATS.NETCDF,
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
                            supported_formats=FORMATS.NETCDF,
                            as_reference=True)

    solution = ComplexOutput('solution', 'solution.rvc file to restart another simulation with the conditions '
                                         'at the end of this simulation.',
                             supported_formats=FORMATS.TEXT,
                             as_reference=True)

    diagnostics = ComplexOutput('diagnostics', 'Performance diagnostic values',
                                abstract="Model diagnostic CSV file.",
                                supported_formats=FORMATS.TEXT,
                                as_reference=True)

    def __init__(self):

        inputs = [self.nc, self.params, self.start_date, self.end_date, self.duration, self.init, self.run_name,
                  self.name, self.area, self.latitude, self.longitude, self.elevation]

        outputs = [self.hydrograph, self.storage, self.solution, self.diagnostics]

        super(RavenGR4JCemaNeigeProcess, self).__init__(
            self._handler,
            identifier=self.identifier,
            title=self.title,
            version=self.version,
            abstract=self.abstract,
            inputs=inputs,
            outputs=outputs,
            status_supported=True,
            store_supported=True
        )

    def _handler(self, request, response):

        # -------------- #
        #  Model config  #
        # -------------- #

        # Default values for this particular process
        rvi, rvp, rvc, rvh, rvt = (self.defaults[k] for k in ['rvi', 'rvp', 'rvc', 'rvh', 'rvt'])

        # Assign the correct input forcing file and variable names fpr each field
        files = [i.file for i in request.inputs['nc']]
        rvt.update(ravenio.assign_files(files, rvt.keys()))

        # Read model configuration and simulation information
        for k, v in rvi.items():
            if v is None and k.lower() in request.inputs.keys():
                rvi[k] = request.inputs[k.lower()][0].data

        # Read basin attributes
        for k, v in rvh.items():
            if v is None:
                rvh[k] = request.inputs[k.lower()][0].data

        # Assemble model configuration parameters
        rvp.update(dict(zip(self.pdefaults.keys(), map(float, request.inputs['params'][0].data.split(',')))))

        # Assemble soil initial conditions
        rvc.update(dict(zip(rvc.keys(), map(float, request.inputs['init'][0].data.split(',')))))

        # Handle start and end date defaults
        start, end = ravenio.start_end_date(files)

        if rvi['Start_Date'] == dt.datetime(1, 1, 1):
            rvi['Start_Date'] = start

        if rvi['Duration'] > 0:
            if rvi['End_Date'] != dt.datetime(1, 1, 1):
                LOGGER.warning("Ambiguous input detected, values for Duration and End_Date have been specified."
                               "Defaults to Duration value.")
        else:
            if rvi['End_Date'] == dt.datetime(1, 1, 1):
                rvi['End_Date'] = end

            rvi['Duration'] = (rvi['End_Date'] - rvi['Start_Date']).days

        # Prepare simulation subdirectory
        params = dict(rvi=rvi, rvp=rvp, rvc=rvc, rvh=rvh, rvt=rvt)
        cmd = ravenio.setup_model(self.identifier, self.workdir, params)

        # Run the simulation
        subprocess.call(cmd)

        # Output files default names
        out_files = {'hydrograph': 'Hydrographs.nc',
                     'storage': 'WatershedStorage.nc',
                     'solution': 'solution.rvc',
                     'diagnostics': 'Diagnostics.csv'}

        # Assign the response outputs to the full names
        sim = rvi['run_name']
        for key, val in out_files.items():
            fn = os.path.join(self.workdir, 'output', '_'.join([sim, val]))
            response.outputs[key].file = fn

        return response
