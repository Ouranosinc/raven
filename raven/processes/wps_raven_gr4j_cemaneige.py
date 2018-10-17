from pywps import Process
from pywps import LiteralInput, LiteralOutput
# from pywps import BoundingBoxInput
from pywps import BoundingBoxOutput
from pywps import ComplexInput, ComplexOutput
from pywps import Format, FORMATS
from pywps.app.Common import Metadata


import logging
LOGGER = logging.getLogger("PYWPS")


defaults = dict(
    rvi=dict(StartDate=None, Duration=None, TimeStep=1.0, EvaluationMetrics='NASH_SUTCLIFFE RMSE'),
    rvp=dict(GR4J_X1=None, GR4J_X2=None, GR4J_X3=None, GR4J_X4=None, AvgAnnualSnow=None, AirSnowCoeff=None),
    rvc=dict(SOIL_0=None, SOIL_1=None),
    rvh=dict(NAME=None, AREA=None, ELEVATION=None, LATITUDE=None, LONGITUDE=None),


class RavenGR4JCemaNeigeProcess(Process):
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
                           min_occurs=0)
              ],
    outputs = [ComplexOutput('q', 'Discharge time series (mm)',
                             supported_formats=[FORMATS.NETCDF],
                             as_reference=True)]

    def __init__(self):
        super(GR4JCemaNeigeProcess, self).__init__(
            self._handler,
            identifier='raven_gr4j_cemaneige',
            title='',
            version='',
            abstract='GR4J + CEMANEIGE hydrological model',
            inputs=inputs,
            outputs=outputs,
            status_supported=True,
            store_supported=True
        )

    def _handler(self, request, response):
        from models import gr4j
        import pandas as pd

        # Read the input data
        pr = xr.open_dataset(request.inputs['pr'][0].url).pr
        tas = xr.open_dataset(request.inputs['tas'][0].url).tas
        evap = xr.open_dataset(request.inputs['evap'][0].url).evap
        sd = request.inputs['start_date'][0].data

        if sd == self.inputs['start_date'].default:
            sd = pr.time[0]

        df = pd.DataFrame({'Temp': tas, 'Evap': evap, 'Prec': pr}).where(df.time >= sd)

        params = map(float, request.inputs['params'][0].data.split(','))

        q = gr4j(df, params)
        da = xr.DataArray(q)
        da.to_netcdf()

