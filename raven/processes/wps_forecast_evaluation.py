import json

from pywps import ComplexInput, ComplexOutput
from pywps import FORMATS
from pywps import LiteralInput
from pywps import Process
from pywps.app.Common import Metadata

import numpy as np
from scipy.stats import norm
import xskillscore as xs
import xarray as xr

class ForecastEvaluationProcess(Process):
    def __init__(self):
        inputs = [ComplexInput('obs', 'Stream flow observation',
                               abstract='Steam flow observation time series',
                               supported_formats=(FORMATS.NETCDF,)),
                  ComplexInput('fcst', 'Stream flow forecast',
                               abstract='Stream flow forecast time series, deterministic or ensemble',
                               supported_formats=(FORMATS.NETCDF,)),
                  LiteralInput('skip_nans', 'Skip NaNs in metric evaluation',
                               abstract='Skip NaNs in forecast evaluation computations',
                               data_type='boolean',
                               default = True,
                               min_occurs=0,
                               max_occurs=1),
                  LiteralInput('BSS_threshold', 'Threshold for Brier Skill Score exceeding given threshold',
                               abstract='Threshold for Brier Skill Score exceeding given threshold',
                               data_type='float',
                               default = 0.7,
                               min_occurs=0,
                               max_occurs=1),                         
                  LiteralInput('name', 'Forecast evaluation metric name',
                               abstract="One or multiple Forecast evaluation metric names. If None, defaults to all.",
                               data_type='string',
                               #allowed_values=tuple(funcs.keys()),
                               default=None,
                               min_occurs=0,
                               max_occurs=17) ### CHANGE THIS TO THE CORRECT NUMBER
                  ]

        outputs = [ComplexOutput('metrics', 'Forecast evaluation metrics values',
                                 abstract="Returns up to ###CHANGE HERE### objective function values, depending on the user's "
                                          "requests. By default all ###CHANGE HERE### are returned. JSON dictionary format.",
                                 supported_formats=(FORMATS.JSON, )),
                   ]

        super(ForecastEvaluationProcess, self).__init__(
            self._handler,
            identifier="forecast-evaluation",
            title="Forecast evaluation based on the XSkillScore package for deterministic and ensemble forecasts.",
            version="1.0",
            abstract="This process takes two NETCDF files (one containing the observed and the other the forecast data) "
                     "and computes forecast evaluation metrics between them. Metrics are calculated according to if there are multiple members in the dataset (Probabilistic) or not (Deterministic)",
            metadata=[Metadata("XSkillScore Documentation", "https://pypi.org/project/xskillscore/")],
            inputs=inputs,
            outputs=outputs,
            #keywords=["forecast evaluation", "ensemble", "deterministic"] + list(funcs.keys()),
            keywords=["forecast evaluation", "ensemble", "deterministic"],
            status_supported=True,
            store_supported=True)

    def _handler(self, request, response):

        
        obs_fn = request.inputs['obs'][0].file
        fcst_fn = request.inputs['fcst'][0].file
        
        skipna=request.inputs['skip_nans'][0].data
        BSS_threshold=request.inputs['BSS_threshold'][0].data
  
        obs = xr.open_dataset(obs_fn)
        fcst = xr.open_dataset(fcst_fn)

        # NaNs are handled by default in XSkillScore
        out = {}
        if not 'member' in fcst:
            obs=obs.to_array()[0]
            fcst=fcst.to_array()[0]

            ### Deterministic metrics
            # Pearson's correlation coefficient
            out['r'] = np.array([xs.pearson_r(obs, fcst, "time", skipna=skipna).data])[0]
            
            # 2-tailed p-value of Pearson's correlation coefficient
            out['r_p_value'] = np.array([xs.pearson_r_p_value(obs, fcst, "time", skipna=skipna).data])[0]
            
            # R^2 (coefficient of determination) score
            out['r2'] = np.array([xs.r2(obs, fcst, "time", skipna=skipna).data])[0]
            
            # Spearman's correlation coefficient
            out['rs'] = np.array([xs.spearman_r(obs, fcst, "time", skipna=skipna).data])[0]
            
            # 2-tailed p-value associated with Spearman's correlation coefficient
            out['rs_p_value'] = np.array([xs.spearman_r_p_value(obs, fcst, "time", skipna=skipna).data])[0]
            
            # Root Mean Squared Error
            out['rmse'] = np.array([xs.rmse(obs, fcst, "time", skipna=skipna).data])[0]
            
            # Mean Squared Error
            out['mse'] = np.array([xs.mse(obs, fcst, "time", skipna=skipna).data])[0]
            
            # Mean Absolute Error
            out['mae'] = np.array([xs.mae(obs, fcst, "time", skipna=skipna).data])[0]
            
            # Median Absolute Error
            out['median_absolute_error'] = np.array([xs.median_absolute_error(obs, fcst, "time", skipna=skipna).data])[0]
            
            # Mean Absolute Percentage Error
            out['mape'] = np.array([xs.mape(obs, fcst, "time", skipna=skipna).data])[0]
            
            # Symmetric Mean Absolute Percentage Error
            out['smape'] = np.array([xs.smape(obs, fcst, "time", skipna=skipna).data])[0]
        
        
        ### Probabilistic
        if 'member' in fcst:
            
            obs=obs.to_array()[0]
            fcst=fcst.to_array()[0]
            # Continuous Ranked Probability Score with the ensemble distribution
            tmp=np.empty(len(obs['time']),dtype=object)
            for i in range(0,len(obs['time'])):
                tmp[i] = xs.crps_ensemble(obs.isel(time=i), fcst.isel(time=i)).data
            out['crps_ensemble']=np.mean(tmp)
     
            # Continuous Ranked Probability Score with a Gaussian distribution
            tmp=np.empty(len(obs['time']),dtype=object)
            for i in range(0,len(obs['time'])):
                tmp[i] = xs.crps_gaussian(obs.isel(time=i), fcst.isel(time=i).mean("member"), fcst.isel(time=i).std("member")).data
            out['crps_gaussian'] = np.mean(tmp)
            
            # Continuous Ranked Probability Score with numerical integration
            # of the normal distribution
            tmp=np.empty(len(obs['time']),dtype=object)
            for i in range(0,len(obs['time'])):
                tmp[i] = xs.crps_quadrature(obs.isel(time=i), norm).data
            out['crps_quadrature'] = np.mean(tmp)
            
            # Brier scores of an ensemble for exceeding given thresholds
            tmp=np.empty(len(obs['time']),dtype=object)
            for i in range(0,len(obs['time'])):
                tmp[i] = xs.threshold_brier_score(obs.isel(time=i), fcst.isel(time=i), BSS_threshold)
            out['threshold_brier_score'] = np.mean(tmp)
            
            # Brier score
            tmp=np.empty(len(obs['time']),dtype=object)
            for i in range(0,len(obs['time'])):
                tmp[i] = xs.brier_score(obs.isel(time=i) > 0.5, (fcst.isel(time=i) > 0.5).mean("member"))
            out['brier_score'] = np.mean(tmp)
           
        response.outputs['metrics'].data = json.dumps(out)
        return response
