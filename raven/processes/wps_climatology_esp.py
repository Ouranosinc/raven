import logging

from pywps import LiteralInput, Process
from ravenpy.utilities import forecasting

from . import wpsio as wio


class ClimatologyEspProcess(Process):
    def __init__(self):

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

        # The number of parameters will depend on the selected model in model_name
        params = LiteralInput(
            "params",
            "Comma separated list of model parameters",
            abstract="Parameters to run the model",
            data_type="string",
            min_occurs=1,
            max_occurs=1,
        )

        inputs = [
            fdate,
            duration,
            params,
            wio.ts,
            wio.latitude,
            wio.longitude,
            wio.name,
            wio.model_name,
            wio.area,
            wio.elevation,
        ]

        outputs = [
            wio.forecast,
        ]

        super(ClimatologyEspProcess, self).__init__(
            self._handler,
            identifier="climatology_esp",
            title="",
            version="1.0",
            abstract="",
            metadata=[],
            inputs=inputs,
            outputs=outputs,
            keywords=[],
            status_supported=True,
            store_supported=True,
        )

    def _handler(self, request, response):

        # Collect info from call
        kwds = {}
        for key, val in request.inputs.items():
            kwds[key] = request.inputs[key][0].data

        kwds["ts"] = request.inputs["ts"][0].file

        # Get info from kwds but remove the ones that are not understood by RAVEN.
        # Have to delete these because or else I get an error: no config key named model_name
        forecast_date = kwds.pop("forecast_date")
        forecast_duration = kwds.pop("forecast_duration")
        model_name = kwds.pop("model_name")

        # Get the model parameters, transform them to a list of floats and write them back to the kwds config.
        params = kwds["params"]
        csv = params.replace("(", "").replace(")", "")
        params = list(map(float, csv.split(",")))
        kwds["params"] = params

        qsims = forecasting.perform_climatology_esp(
            model_name, forecast_date, forecast_duration, workdir=self.workdir, **kwds
        )

        # Prepare the forecast netcdf result file and send the path to the results output.
        forecastfile = self.workdir + "/forecast.nc"
        qsims.to_netcdf(forecastfile)
        response.outputs["forecast"].file = forecastfile

        return response
