from pywps import ComplexInput, ComplexOutput, LiteralInput
from pywps import FORMATS, Format
from pywps import Process

import xarray as xr
import tempfile


class ForecastFloodRiskProcess(Process):
    def __init__(self):
        inputs = [
            ComplexInput(
                "fcst",
                "Streamflow forecasts (ensemble or deterministic)",
                abstract="NetCDF file including streamflow forecasts (ensemble or deterministic) to evaluate the flood risk based on the flood level (threshold).",
                supported_formats=(FORMATS.NETCDF,),
            ),
            LiteralInput(
                "name",
                "Name of the variable within the fcst netcdf file.",
                abstract="Name of the variable within the fcst netcdf file to access the data using a variable naming scheme",
                data_type="string",
                min_occurs=1,
                max_occurs=1,
            ),
            LiteralInput(
                "flood_level",
                "Flood level threshold.",
                abstract="Flood level threshold. Will be used to determine if forecasts exceed this specified flood threshold.",
                data_type="float",
                min_occurs=1,
                max_occurs=1,
            ),
        ]

        outputs = [
            ComplexOutput(
                "flood_risk",
                "Empirical flood risk based on forecast data and specified flood level.",
                supported_formats=[
                    FORMATS.NETCDF,
                    Format("application/zip", extension=".zip", encoding="base64"),
                ],
                abstract="Netcdf file including ratio of members that exceed the flood level for each forecast day",
                as_reference=True,
            )
        ]

        super(ForecastFloodRiskProcess, self).__init__(
            self._handler,
            identifier="forecast-floodrisk",
            title="Calculate flood risk from a forecast and flood level threshold",
            version="1.0",
            abstract="This WPS service takes a NetCDF forecast (ensemble or deterministic) and returns the empirical exceedance probability for each forecast day based on a flood level threshold. "
            "Results are returned in a NetCDF file with a time series of probabilities of flood level exceedance.",
            inputs=inputs,
            outputs=outputs,
            keywords=[
                "forecast",
                "flood",
                "risk",
                "threshold",
                "ensemble",
                "deterministic",
            ],
            status_supported=True,
            store_supported=True,
        )

    def _handler(self, request, response):

        # Read inputs from request
        fcst_fn = request.inputs["fcst"][0].file
        fcst_var = request.inputs["name"][0].data
        flood_level = request.inputs["flood_level"][0].data

        # Open netCDF files
        fcst_ds = xr.open_dataset(fcst_fn)

        # Get variable names
        fcst = fcst_ds[fcst_var]

        # ---- Calculations ---- #
        # Ensemble: for each day, calculate the percentage of members that are above the threshold
        if "member" in fcst.coords:
            # Get number of members originally
            number_members = len(fcst.member)

            # now compute the ratio of cases that are above the threshold
            pct = (
                fcst.where(fcst > flood_level)
                .notnull()
                .sum(dim="member", keep_attrs=True)
                / number_members
            )

        # it's deterministic:
        else:
            pct = (
                fcst.where(fcst > flood_level).notnull() / 1.0
            )  # This is needed to return values instead of floats

        # Put in file
        filename = tempfile.mkdtemp() + "/fcst_flood_risk.nc"
        pct.to_netcdf(filename)

        # Return the response.
        response.outputs["flood_risk"].file = str(filename)
        return response
