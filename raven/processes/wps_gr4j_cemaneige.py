from pywps import Process
from pywps import LiteralInput
from pywps import ComplexInput, ComplexOutput
from pywps import FORMATS
from pywps.app.Common import Metadata
import xarray as xr
import os
import datetime as dt

import logging
LOGGER = logging.getLogger("PYWPS")


class GR4JCemaNeigeProcess(Process):
    def __init__(self):
        inputs = [ComplexInput('pr', 'Precipitation (mm)',
                               abstract='netCDF file storing daily precipitation time series (pr).',
                               min_occurs=1,
                               supported_formats=[FORMATS.NETCDF]),

                  ComplexInput('tas', 'Temperature (K)',
                               abstract='netCDF file storing daily temperature time series (tas).',
                               min_occurs=1,
                               supported_formats=[FORMATS.NETCDF]),

                  ComplexInput('evap', 'Evapotranspiration (mm)',
                               abstract='netCDF file storing daily evapotranspiration time series (evap).',
                               min_occurs=1,
                               supported_formats=[FORMATS.NETCDF]),

                  LiteralInput('params', 'Comma separated list of model parameters',
                               abstract='Parameters: GR4J_X1, GR4J_X2, GR4J_X3, GR4J_X4, CN_X1, 1-CN_X2'
                                        'Raven: SOIL_PROD, GR4J_X2, GR4J_X3, GR4J_X4, AvgAnnualSnow, AirSnowCoeff',
                               min_occurs=0,
                               data_type='string',
                               default='300, -1, 200, 2, .5, 8'),

                  LiteralInput('start_date', 'Simulation start date (AAAA-MM-DD)',
                               abstract='Start date of the simulation (AAAA-MM-DD). '
                                        'Defaults to the start of the forcing file. ',
                               data_type='date',
                               default='0001-01-01',
                               min_occurs=0),

                  LiteralInput('end_date', 'Simulation end date (AAAA-MM-DD)',
                               abstract='End date of the simulation (AAAA-MM-DD). '
                                        'Defaults to the end of the forcing file.',
                               data_type='date',
                               default='0001-01-01',
                               min_occurs=0),
                  ]

        outputs = [ComplexOutput('q', 'Discharge time series (mm)',
                                 supported_formats=[FORMATS.NETCDF],
                                 as_reference=True), ]

        super(GR4JCemaNeigeProcess, self).__init__(
            self._handler,
            identifier='gr4j_cemaneige',
            title='',
            version='',
            abstract='GR4J + CEMANEIGE hydrological model',
            inputs=inputs,
            outputs=outputs,
            status_supported=True,
            store_supported=True
        )

    def _handler(self, request, response):
        from raven.models import gr4j

        xr.set_options(enable_cftimeindex=True)

        # Read the input data
        paths = [request.inputs['pr'][0].file,
                 request.inputs['tas'][0].file,
                 request.inputs['evap'][0].file]
        ds = xr.open_mfdataset(paths).rename({'pr': 'Prec', 'tas': 'Temp', 'evap': 'Evap'})

        # Convert temperature to Celsius
        ds['Temp'] -= 273.15

        # Convert to pandas.DataFrame
        df = ds.to_dataframe()

        # Parse start and end time
        start = request.inputs['start_date'][0].data
        end = request.inputs['end_date'][0].data

        if start == dt.date(1, 1, 1):
            start = df.index[0]
        else:  # Convert date to datetime
            start = dt.datetime.combine(start, dt.time())

        if end == dt.date(1, 1, 1):
            end = df.index[-1]
        else:
            end = dt.datetime.combine(end, dt.time())

        # Select time steps
        mask = (df.index >= start) & (df.index <= end)
        df = df.loc[mask]

        # Create the model parameter array
        params = map(float, request.inputs['params'][0].data.split(','))

        # Run the model
        q = gr4j(df, params)

        # Set output attributes
        da = xr.DataArray(q)
        da.attrs['long_name'] = 'discharge'
        da.attrs['standard_name'] = 'discharge'
        da.attrs['units'] = 'mm/day'

        # Save to disk
        fn = os.path.join(self.workdir, 'q.nc')
        da.to_netcdf(fn)

        # Store file path in process response
        response.outputs['q'].file = fn

        return response
