from pywps import Process
from pywps import LiteralInput, LiteralOutput
# from pywps import BoundingBoxInput
from pywps import BoundingBoxOutput
from pywps import ComplexInput, ComplexOutput
from pywps import Format, FORMATS
from pywps.app.Common import Metadata
import xarray as xr
import os


import logging
LOGGER = logging.getLogger("PYWPS")


class GR4JCemaNeigeProcess(Process):
    def __init__(self):
        inputs = [ComplexInput('pr', 'Precipitation (mm)',
                               abstract='netCDF file storing daily precipitation time series (pr).',
                               min_occurs=1,
                               supported_formats=[FORMATS.NETCDF, FORMATS.DODS]),
                  ComplexInput('tas', 'Temperature (K)',
                               abstract='netCDF file storing daily temperature time series (tas).',
                               min_occurs=1,
                               supported_formats=[FORMATS.NETCDF, FORMATS.DODS]),
                  ComplexInput('evap', 'Evapotranspiration (mm)',
                               abstract='netCDF file storing daily evapotranspiration time series (evap).',
                               min_occurs=1,
                               supported_formats=[FORMATS.NETCDF, FORMATS.DODS]),
                  LiteralInput('params', 'Comma separated list of model parameters',
                               abstract='Parameters: GR4J_X1, GR4J_X2, GR4J_X3, GR4J_X4, CN_X1, 1-CN_X2'
                                        'Raven: SOIL_PROD, GR4J_X2, GR4J_X3, GR4J_X4, AvgAnnualSnow, AirSnowCoeff',
                               min_occurs=0,
                               default='0.696, 0.7, 19.7, 2.09, 123.3, 0.75'),
                  LiteralInput('start_date', 'Simulation start date',
                               abstract='Start date of the simulation. Defaults to the start of the forcing file.',
                               data_type='datetime',
                               default='01-01-0001',
                               min_occurs=0),
                  LiteralInput('end_date', 'Simulation end date',
                               abstract='End date of the simulation. Defaults to the end of the forcing file.',
                               data_type='datetime',
                               default='01-01-0001',
                               min_occurs=0)
                  ],
        outputs = [ComplexOutput('q', 'Discharge time series (mm)',
                                 supported_formats=[FORMATS.NETCDF],
                                 as_reference=True)]

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
        import pandas as pd

        # Read the input data
        pr = xr.open_dataset(request.inputs['pr'][0].url).pr
        tas = xr.open_dataset(request.inputs['tas'][0].url).tas
        evap = xr.open_dataset(request.inputs['evap'][0].url).evap

        start = request.inputs['start_date'][0].data
        end = request.inputs['end_date'][0].data

        if start == self.inputs['start_date'].default:
            start = pr.time[0]

        if end == self.inputs['end_date'].default:
            end = pr.time[-1]

        # Create the input dataframe
        df = pd.DataFrame({'Temp': tas, 'Evap': evap, 'Prec': pr})
        df = df.where(df.time >= start and df.time < end)  # really ?

        # Create the model parameter array
        params = map(float, request.inputs['params'][0].data.split(','))

        # Run the model
        q = gr4j(df, params)

        da = xr.DataArray(q)
        da.attrs['long_name'] = 'discharge'
        da.attrs['standard_name'] = 'discharge'
        da.attrs['units'] = 'mm/day'

        fn = os.path.join(self.workdir, 'q.nc')
        da.to_netcdf(fn)
        response.outputs['q'].file = fn
