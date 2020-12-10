import logging
from pathlib import Path
import fiona
import xarray as xr
import pandas as pd
from raven.utilities import forecasting
from raven.models import RavenMultiModel
from . wps_raven_multi_model import RavenMultiModelProcess, hmets, gr4jcn, hbvec
from . import wpsio as wio

LOGGER = logging.getLogger("PYWPS")


class RealtimeForecastProcess(RavenMultiModelProcess):
    identifier = "realtime-forecast"
    title = "Perform realtime forecast using most recent ECCC forecast data."
    abstract = "Perform a deterministic or probabilistic raven forecast using the most recent ECCC forecast."
    version = '0.1'
    keywords = ["forecasting", "ECCC", "GEPS", "REPS", "GDPS", "RDPS", "ensemble forecasts"]

    inputs = [wio.forecast_model, wio.region_vector, hmets, gr4jcn, hbvec, wio.duration,
              wio.run_name, wio.name, wio.area, wio.latitude, wio.longitude, wio.elevation, wio.rain_snow_fraction,
              wio.nc_spec, wio.rvc]
    model_cls = RavenMultiModel

    def meteo(self, request):
        """Fetch the latest forecast from ECCC.

        Returns
        -------
        list
          List of input file paths.
        """
        # Prepare the forecast data from the latest ECCC forecast.

        # Use for testing
        #return [Path("~/src/raven-testdata/eccc_geps/fcstfile.nc"),]

        # Region shapefile
        vector_file = self.region(request)

        # Forecast model
        forecast_model = request.inputs.pop("forecast_model")[0].data

        # Fetch data and average over region
        fcst = forecasting.get_recent_ECCC_forecast(
               fiona.open(vector_file), climate_model=forecast_model)

        # Write the forecast data to file on-disk
        # To speed-up testing, copy this file and return it instead of recomputing every time.
        fn = Path(self.workdir) / "fcstfile.nc"
        fcst.to_netcdf(fn)

        return [fn]

    def run(self, model, ts, kwds):
        """Initialize the model with the RVC file, then run it with the forecast data."""
        # Open forecast file and set some attributes.
        fcst = xr.open_dataset(ts[0])
        kwds['nc_index'] = range(fcst.dims.get("member"))

        model(ts=ts, **kwds)
        fcst.close()


