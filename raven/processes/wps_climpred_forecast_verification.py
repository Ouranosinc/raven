from pathlib import Path

import xarray as xr
from pywps import FORMATS, ComplexInput, ComplexOutput, LiteralInput, Process
from ravenpy.utilities.forecasting import make_climpred_hindcast_object


class ClimpredHindcastVerificationProcess(Process):
    def __init__(self):

        hindcasts = ComplexInput(
            "hindcasts",
            "3-dimensional xarray dataset / netcdf with hindcasts",
            abstract="The 3D netcdf dataset that contains the init, member and lead variables",
            supported_formats=[FORMATS.NETCDF],
            min_occurs=1,
            max_occurs=1,
        )

        observations = ComplexInput(
            "observations",
            "1-dimensional xarray dataset / netcdf with flow observations",
            abstract="The 1D netcdf with the observed streamflow for verification",
            supported_formats=[FORMATS.NETCDF],
            min_occurs=1,
            max_occurs=1,
        )

        metric = LiteralInput(
            "metric",
            'Verification metric. Can be ["rank_histogram","crps" or "reliability"]',
            data_type="string",
            abstract='Name of the verification metric. Can be ["rank_histogram","crps" or "reliability"]',
            min_occurs=1,
            max_occurs=1,
        )

        inputs = [hindcasts, observations, metric]

        outputs = [
            ComplexOutput(
                "verification_metrics",
                "The verification_metrics dataset as computed by climpred, ready to plot.",
                supported_formats=[FORMATS.NETCDF],
                abstract="Netcdf file including the verification metrics that can be used for plotting hindcast performance. Contents vary according to input metric",
                as_reference=True,
            )
        ]

        super(ClimpredHindcastVerificationProcess, self).__init__(
            self._handler,
            identifier="climpred_hindcast_verification",
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

        hindcasts = xr.open_dataset(request.inputs["hindcasts"][0].file)
        qobs = xr.open_dataset(request.inputs["observations"][0].file)
        metric = request.inputs["metric"][0].data

        # Here units are days and always will be!
        hindcasts["lead"] = list(range(1, len(hindcasts.lead) + 1))
        hindcasts["lead"].attrs[
            "units"
        ] = "days"  # Note that this should be done in the make_climpred function but does not follow after to_netcdf.

        # Once we have the correctly formatted datasets, Make the hindcast object for climpred
        hindcast_object = make_climpred_hindcast_object(hindcasts, qobs)

        # This function is used to convert to binary to see if yes/no forecast is larger than obs
        def pos(x):
            return x > 0  # Check for binary outcome

        # These three functions respectively compute the rank histogram,
        # the crps and the reliability for the set of initialized dates
        # (i.e. forecast issue dates, here 1 day per year at the same calendar day).
        if metric == "rank_histogram":
            verification_metrics = hindcast_object.verify(
                metric="rank_histogram",
                comparison="m2o",
                dim=["member", "init"],
                alignment="same_inits",
            )
        if metric == "crps":
            verification_metrics = hindcast_object.verify(
                metric="crps",
                comparison="m2o",
                dim=["member", "init"],
                alignment="same_inits",
            )
        if metric == "reliability":
            verification_metrics = hindcast_object.verify(
                metric="reliability",
                comparison="m2o",
                dim=["member", "init"],
                alignment="same_inits",
                logical=pos,
            )

        assert "flow" in verification_metrics

        assert verification_metrics.flow.shape[0] == hindcasts.lead.shape[0]

        verif_metrics_file = Path(self.workdir) / "verification_metrics.nc"
        verification_metrics.to_netcdf(verif_metrics_file)
        response.outputs["verification_metrics"].file = str(verif_metrics_file)

        return response
