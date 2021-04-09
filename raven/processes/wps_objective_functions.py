import json

import spotpy as sp
import xarray as xr
from pywps import FORMATS, ComplexInput, ComplexOutput, LiteralInput, Process
from pywps.app.Common import Metadata

funcs = {f.__name__: f for f in sp.objectivefunctions._all_functions}


class ObjectiveFunctionProcess(Process):
    def __init__(self):
        inputs = [
            ComplexInput(
                "obs",
                "Stream flow observation",
                abstract="Steam flow observation time series",
                supported_formats=(FORMATS.NETCDF,),
            ),
            ComplexInput(
                "sim",
                "Stream flow simulation",
                abstract="Stream flow simulation time series",
                supported_formats=(FORMATS.NETCDF,),
            ),
            LiteralInput(
                "name",
                "Objective function name",
                abstract="One or multiple objective function name. If None, defaults to all.",
                data_type="string",
                allowed_values=tuple(funcs.keys()),
                default=None,
                min_occurs=0,
                max_occurs=17,
            ),
        ]

        outputs = [
            ComplexOutput(
                "metrics",
                "Objective function values",
                abstract="Returns up to 17 objective function values, depending on the user's "
                "requests. By default all 17 are returned. JSON dictionary format.",
                supported_formats=(FORMATS.JSON,),
            ),
        ]

        super(ObjectiveFunctionProcess, self).__init__(
            self._handler,
            identifier="objective-function",
            title="Objective-function process based on SpotPy and its 17 objective functions.",
            version="1.0",
            abstract="This process takes two NETCDF files (one containing variable 'q_sim' and the other 'q_obs') "
            "and computes objective-function metrics between them.",
            metadata=[
                Metadata(
                    "SPOTPY Documentation",
                    "http://fb09-pasig.umwelt.uni-giessen.de/spotpy/",
                )
            ],
            inputs=inputs,
            outputs=outputs,
            keywords=["objective functions", "hydrological signatures", "optimization"]
            + list(funcs.keys()),
            status_supported=True,
            store_supported=True,
        )

    def _handler(self, request, response):

        obs_fn = request.inputs["obs"][0].file
        sim_fn = request.inputs["sim"][0].file
        if "name" in request.inputs:
            names = [i.data for i in request.inputs["name"]]
        else:
            names = funcs.keys()

        obs = xr.open_dataset(obs_fn)
        sim = xr.open_dataset(sim_fn)

        # There is no support yet for handling NaN in SpotPy. Here we're starting from the second index to avoid missing
        # values in the first time index for obs.
        out = {}
        for name in names:
            da = xr.apply_ufunc(
                funcs[name],
                obs["q_obs"].isel(time=slice(1, None)),
                sim["q_sim"].isel(time=slice(1, None)),
                input_core_dims=[
                    [
                        "time",
                    ],
                    [
                        "time",
                    ],
                ],
                vectorize=True,
            )

            # For now we're assuming there is just one basin
            out[name] = da.data[0]

        response.outputs["metrics"].data = json.dumps(out)
        return response
