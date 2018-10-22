import datetime as dt
import os
from pywps import Process
from pywps import LiteralInput, LiteralOutput
from pywps import ComplexInput, ComplexOutput
from pywps import Format, FORMATS
from pywps.app.Common import Metadata
import subprocess
from . import ravenio

import logging
from collections import OrderedDict as Odict

LOGGER = logging.getLogger("PYWPS")


"""
Notes
-----

The configuration files for RAVEN's GR4J-Cemaneige model and in models/raven-gr4j. 
All parameters that could potentially be user-defined are tagged using {}.
Different WPS processes can provide different level of customization for the same model. 
The idea is to use the `defaults` dictionary to set *frozen* parameters.   
The Process itself can also set defaults for convenience. 
"""


defaults = Odict(
    rvi=dict(Start_Date=None, End_Date=None, Duration=None, TimeStep=1.0, EvaluationMetrics='NASH_SUTCLIFFE RMSE'),
    rvp=Odict(GR4J_X1=None, GR4J_X2=None, GR4J_X3=None, GR4J_X4=None, AvgAnnualSnow=None, AirSnowCoeff=None),
    rvc=Odict(SOIL_0=None, SOIL_1=None),
    rvh=dict(NAME=None, AREA=None, ELEVATION=None, LATITUDE=None, LONGITUDE=None),
    rvt=dict(RAIN=None, SNOW=None, TMIN=None, TMAX=None, PET=None, QOBS=None)
)

class RavenGR4JCemaNeigeProcess(Process):

    def __init__(self):

        inputs = [ComplexInput('nc', 'netCDF input files',
                               abstract='NetCDF file or files storing'
                                        ' daily liquid precipitation (rain [mm]), '
                                        'solid precipitation (snow [mm]), '
                                        'minimum temperature (tasmin [degC]), '
                                        'maximum temperature (tasmax [degC]), '
                                        'potential evapotranspiration (pet [mm]) and '
                                        'observed streamflow (qobs [m3/s]).',
                               min_occurs=1,
                               supported_formats=[FORMATS.NETCDF]),

                  LiteralInput('params', 'Comma separated list of model parameters',
                               abstract='Parameters: SOIL_PROD, GR4J_X2, GR4J_X3, GR4J_X4, AvgAnnualSnow, AirSnowCoeff',
                               data_type='string',
                               default='0.696, 0.7, 19.7, 2.09, 123.3, 0.75',
                               min_occurs=0),

                  LiteralInput('start_date', 'Simulation start date (AAAA-MM-DD)',
                               abstract='Start date of the simulation (AAAA-MM-DD). '
                                        'Defaults to the start of the forcing file. ',
                               data_type='dateTime',
                               default='0001-01-01 00:00:00',
                               min_occurs=0),

                  LiteralInput('end_date', 'Simulation end date (AAAA-MM-DD)',
                               abstract='End date of the simulation (AAAA-MM-DD). '
                                        'Defaults to the end of the forcing file.',
                               data_type='dateTime',
                               default='0001-01-01 00:00:00',
                               min_occurs=0),

                  LiteralInput('duration', 'Simulation duration (days)',
                               abstract='Number of simulated days, defaults to the length of the input forcings.',
                               data_type='nonNegativeInteger',
                               default=0,
                               min_occurs=0),

                  LiteralInput('init', 'Initial soil conditions',
                               abstract='Underground reservoir levels: SOIL_0, SOIL_1',
                               data_type='string',
                               default='0, 0',
                               min_occurs=0),

                  LiteralInput('name', 'Simulation name',
                               abstract='The name given to the simulation, for example <watershed>_<experiment>',
                               data_type='string',
                               default='raven-gr4j-cemaneige-sim',
                               min_occurs=0),

                  LiteralInput('area', 'Watershed area (km2)',
                               abstract='Watershed area (km2)',
                               data_type='float',
                               default=0.,
                               min_occurs=0),

                  LiteralInput('latitude', 'Latitude',
                               abstract="Watershed's centroid latitude",
                               data_type='float',
                               min_occurs=1),

                  LiteralInput('longitude', 'Longitude',
                               abstract="Watershed's centroid longitude",
                               data_type='float',
                               min_occurs=1),

                  LiteralInput('elevation', 'Elevation (m)',
                               abstract="Watershed's mean elevation (m)",
                               data_type='float',
                               min_occurs=1),

                  ]

        outputs = [ComplexOutput('q', 'Discharge time series (mm)',
                                 supported_formats=[FORMATS.NETCDF],
                                 as_reference=True), ]

        super(RavenGR4JCemaNeigeProcess, self).__init__(
            self._handler,
            identifier='raven-gr4j-cemaneige',
            title='',
            version='',
            abstract='Raven GR4J + CEMANEIGE hydrological model',
            inputs=inputs,
            outputs=outputs,
            status_supported=True,
            store_supported=True
            )

    def _handler(self, request, response):

        # -------------- #
        #  Model config  #
        # -------------- #
        rvi, rvp, rvc, rvh, rvt = defaults.values()

        # Assign the correct input forcing file to each field
        files = [i.file for i in request.inputs['nc']]
        vals = ravenio.assign_files(files, [k.lower() for k in rvt.keys()])
        rvt = dict(zip(rvt.keys(), vals))

        # Read model configuration and simulation information
        for k, v in rvi.items():
            if v is None and k.lower() in request.inputs.keys():
                rvi[k] = request.inputs[k.lower()][0].data

        # Read basin attributes
        for k, v in rvh.items():
            if v is None:
                rvh[k] = request.inputs[k.lower()][0].data

        # Assemble model configuration parameters
        rvp.update(dict(zip(rvp.keys(), map(float, request.inputs['params'][0].data.split(',')))))

        # Assemble soil initial conditions
        rvc.update(dict(zip(rvc.keys(), map(float, request.inputs['init'][0].data.split(',')))))
        # -------------- #

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
        subprocess.run(cmd)
        #response.outputs['q'].file = os.path.join(self.workdir, 'output', 'duh.nc')

        return response
