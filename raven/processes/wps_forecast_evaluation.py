import json

from pywps import ComplexInput, ComplexOutput
from pywps import FORMATS
from pywps import LiteralInput
from pywps import Process
from pywps.app.Common import Metadata

import xskillscore as xs
import xarray as xr

# Name of all available metrics
all_metrics = xs.core.deterministic.__all__ + xs.core.probabilistic.__all__

# Remove unsupported metrics
all_metrics.remove('crps_quadrature')
all_metrics.remove("crps_gaussian")
all_metrics.remove("brier_score")


# TODO: Report multidimensional bug for pearson_r to xskillscore.
class HindcastEvaluationProcess(Process):
    def __init__(self):
        inputs = [ComplexInput('obs', 'Stream flow observation',
                               abstract='Stream flow observation time series.',
                               supported_formats=(FORMATS.NETCDF,)),
                  LiteralInput("obs_var", "Observation variable name",
                               abstract="Name of the variable in the observation dataset.",
                               data_type="string",
                               default="q_obs",
                               min_occurs=0,
                               max_occurs=1),
                  ComplexInput('hcst', 'Stream flow hindcast',
                               abstract='Stream flow hindcast time series, deterministic or ensemble.',
                               supported_formats=(FORMATS.NETCDF,)),
                  LiteralInput("hcst_var", "Hindcast variable name",
                               abstract="Name of the variable in the hindcast dataset.",
                               data_type="string",
                               default="q_sim",
                               min_occurs=0,
                               max_occurs=1),
                  LiteralInput('skip_nans', 'Skip NaNs in metric evaluation',
                               abstract='Skip NaNs in hindcast evaluation computations',
                               data_type='boolean',
                               default=True,
                               min_occurs=0,
                               max_occurs=1),
                  LiteralInput('BSS_threshold', 'Threshold for Brier Skill Score exceeding given threshold',
                               abstract='Threshold for Brier Skill Score exceeding given threshold',
                               data_type='float',
                               default=0.7,
                               min_occurs=0,
                               max_occurs=1),
                  LiteralInput('metric', 'Forecast evaluation metric name',
                               abstract="One or multiple hindcast evaluation metric names. If None, defaults to all.",
                               data_type='string',
                               allowed_values=all_metrics,
                               default=None,
                               min_occurs=0,
                               max_occurs=len(all_metrics))
                  ]

        outputs = [ComplexOutput('metrics', 'Hindcast evaluation metrics values',
                                 abstract="JSON dictionary of evaluation metrics averaged over the full period and "
                                          "all members.",
                                 supported_formats=(FORMATS.JSON, )),
                   ]

        super(HindcastEvaluationProcess, self).__init__(
            self._handler,
            identifier="hindcast-evaluation",
            title="Hindcast evaluation based on the XSkillScore package for deterministic and ensemble hindcasts.",
            version="1.0",
            abstract="This process takes two NETCDF files (one containing the observed and the other the hindcast "
                     "data) "
                     "and computes hindcast evaluation metrics between them. Metrics are calculated according to if "
                     "there are multiple members in the dataset (probabilistic) or not (deterministic)",
            metadata=[Metadata("XSkillScore Documentation", "https://pypi.org/project/xskillscore/")],
            inputs=inputs,
            outputs=outputs,
            keywords=["forecast evaluation", "ensemble", "deterministic"] + list(all_metrics),
            status_supported=True,
            store_supported=True)

    def _handler(self, request, response):

        # Read inputs from request
        obs_fn = request.inputs['obs'][0].file
        obs_var = request.inputs['obs_var'][0].data

        hcst_fn = request.inputs['hcst'][0].file
        hcst_var = request.inputs['hcst_var'][0].data

        skipna = request.inputs['skip_nans'][0].data
        bss_threshold = request.inputs['BSS_threshold'][0].data

        if 'metric' in request.inputs:
            metrics = [m.data for m in request.inputs['metric']]
        else:
            metrics = all_metrics

        # Open netCDF files
        obs_ds = xr.open_dataset(obs_fn)
        hcst_ds = xr.open_dataset(hcst_fn)

        # Get variable names
        obs = obs_ds[obs_var]
        hcst = hcst_ds[hcst_var]

        # ---- Calculations ---- #
        # NaNs are handled by default in XSkillScore
        out = {}
        for metric in metrics:
            func = getattr(xs, metric)

            if metric in xs.core.deterministic.__all__:
                m = func(obs, hcst, dim="time", skipna=skipna)

            elif "member" in hcst.dims:

                if metric == "threshold_brier_score":
                    m = func(obs, hcst, threshold=bss_threshold).mean("time")

                else:
                    m = func(obs, hcst, dim="member").mean("time")

            out[metric] = m.values.tolist()

        response.outputs['metrics'].data = json.dumps(out)
        return response
