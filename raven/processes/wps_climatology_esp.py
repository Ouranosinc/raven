from pywps import LiteralInput
from raven.processes import RavenProcess
from ravenpy.utilities import forecasting
from ravenpy.models import HMETS, GR4JCN, HBVEC, MOHYSE
from . import wpsio as wio

fdate = LiteralInput(
    "forecast_date",
    "date of the forecast",
    abstract="Date that the climatology-based ESP ensemble will be performed",
    data_type="dateTime",
    min_occurs=1,
    max_occurs=1,
)

duration = LiteralInput(
    "forecast_duration",
    "Forecast duration",
    abstract="Duration of the forecast in days",
    data_type="integer",
    default=30,
    min_occurs=0,
    max_occurs=1,
)

params = LiteralInput(
    "params",
    "Comma separated list of model parameters",
    abstract="Parameters to run the model",
    data_type="string",
    min_occurs=1,
    max_occurs=1,
)


class ClimatologyEspProcess(RavenProcess):
    identifier = "climatology_esp"
    abstract = "Forecast based on observed climatology"
    title = "Climatological Ensemble Streamflow Prediction"
    version = "0.1"

    inputs = [
        fdate,
        duration,
        wio.ts,
        wio.model_name,
        params,
        wio.name,
        wio.area,
        wio.latitude,
        wio.longitude,
        wio.elevation,
        wio.nc_spec,
        wio.rvc,
        wio.rain_snow_fraction,
    ]

    outputs = [wio.forecast]
    tuple_inputs = {
        "hmets": HMETS.params,
        "gr4jcn": GR4JCN.params,
        "hbvec": HBVEC.params,
        "mohyse": MOHYSE.params,
    }

    def parse_tuple(self, obj, model_name):
        csv = obj.data.replace("(", "").replace(")", "")
        arr = map(float, csv.split(","))
        return self.tuple_inputs[model_name.lower()](*arr)

    def _handler(self, request, response):
        response.update_status("PyWPS process {} started.".format(self.identifier), 0)

        # Extract params to skip default processing in `self.options`
        params = request.inputs.pop("params")[0]

        # Input data files.
        ts = self.meteo(request)

        # Initial state
        if "rvc" in request.inputs:
            rvc = request.inputs.pop("rvc")[0].file
        else:
            solution = self.get_config(request, ids=("rvc",))
            rvc = list(solution.values()).pop()["rvc"]

        # Model options
        kwds = self.options(request)
        kwds["params"] = self.parse_tuple(params, kwds["model_name"].lower())
        kwds["rvc"] = rvc
        kwds["workdir"] = self.workdir

        # Make forecasts. Only support one input file for now.
        qsims = forecasting.perform_climatology_esp(ts=ts[0], **kwds)

        # Prepare the forecast netcdf result file and send the path to the results output.
        forecastfile = self.workdir + "/forecast.nc"
        qsims.to_netcdf(forecastfile)
        response.outputs["forecast"].file = forecastfile

        return response
