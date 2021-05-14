import logging
from pathlib import Path

import fiona
import xarray as xr
from ravenpy.models import GR4JCN, HBVEC, HMETS, MOHYSE, get_model
from ravenpy.utilities import forecasting

from . import wpsio as wio
from .wps_raven import RavenProcess

LOGGER = logging.getLogger("PYWPS")


class RealtimeForecastProcess(RavenProcess):
    identifier = "realtime-forecast"
    title = "Perform realtime forecast using most recent ECCC forecast data."
    abstract = "Perform a deterministic or probabilistic raven forecast using the most recent ECCC forecast."
    version = "0.1"
    keywords = [
        "forecasting",
        "ECCC",
        "GEPS",
        "REPS",
        "GDPS",
        "RDPS",
        "ensemble forecasts",
    ]
    tuple_inputs = {
        "hmets": HMETS.Params,
        "gr4jcn": GR4JCN.Params,
        "hbvec": HBVEC.Params,
        "mohyse": MOHYSE.Params,
    }
    inputs = [
        wio.forecast_model,
        wio.region_vector,
        wio.hmets,
        wio.gr4jcn,
        wio.hbvec,
        wio.duration,
        wio.run_name,
        wio.area,
        wio.latitude,
        wio.longitude,
        wio.elevation,
        wio.rain_snow_fraction,
        wio.nc_spec,
        wio.rvc,
    ]

    def model(self, request):
        """Return model class."""
        models = list(set(request.inputs.keys()).intersection(self.tuple_inputs.keys()))
        if len(models) > 1:
            raise NotImplementedError("Multi-model simulations are not supported. ")
        name = models.pop()
        params = self.parse_tuple(request.inputs.pop(name)[0])
        model = get_model(name)(workdir=self.workdir)
        model.config.update("params", params)
        return model

    def meteo(self, request):
        """Fetch the latest forecast from ECCC.

        Returns
        -------
        list
          List of input file paths.
        """
        # Region shapefile from request
        vector_file = self.region(request)

        # Forecast model from request
        forecast_model = request.inputs.pop("forecast_model")[0].data

        # Short-cut for testing
        # return [Path("/home/david/src/raven-testdata/eccc_geps/fcstfile.nc"),]

        # Fetch data and average over region
        fcst = forecasting.get_recent_ECCC_forecast(
            fiona.open(vector_file), climate_model=forecast_model
        )

        # Write the forecast data to file on-disk
        # To speed-up testing, copy this file and return it instead of recomputing every time.
        fn = Path(self.workdir) / "fcstfile.nc"
        fcst.to_netcdf(fn)

        return [fn]

    def run(self, model, ts, kwds):
        """Initialize the model with the RVC file, then run it with the forecast data."""
        # Open forecast file and set some attributes.
        fcst = xr.open_dataset(ts[0])
        kwds["nc_index"] = range(fcst.dims.get("member"))

        model(ts=ts, **kwds)

        fcst.close()
