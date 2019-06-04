import json
import logging
from pathlib import Path

from pywps import ComplexInput, FORMATS
from pywps import LiteralInput

from raven.utilities import regionalize, read_gauged_properties, read_gauged_params
from . import wpsio as wio
from .wps_raven import RavenProcess

LOGGER = logging.getLogger("PYWPS")


class RegionalisationProcess(RavenProcess):
    """
    Notes
    -----
    The available regionalization methods are:

    .. glossary::

        Multiple linear regression (MLR)
            Ungauged catchment parameters are estimated individually by a linear regression
            against catchment properties.

        Spatial proximity (SP)
            The ungauged hydrograph is an average of the `n` closest catchments' hydrographs.

        Physical similarity (PS)
            The ungauged hydrograph is an average of the `n` most similar catchments' hydrographs.

        Spatial proximity with inverse distance weighting (SP_IDW)
            The ungauged hydrograph is an average of the `n` closest catchments' hydrographs, but
            the average is weighted using inverse distance weighting

        Physical similarity with inverse distance weighting (PS_IDW)
            The ungauged hydrograph is an average of the `n` most similar catchments' hydrographs, but
            the average is weighted using inverse distance weighting

        Spatial proximity with IDW and regression-based augmentation (SP_IDW_RA)
            The ungauged hydrograph is an average of the `n` closest catchments' hydrographs, but
            the average is weighted using inverse distance weighting. Furthermore, the method uses the CANOPEX/USGS
            dataset to estimate model parameters using Multiple Linear Regression. Parameters whose regression r-squared
            is higher than 0.5 are replaced by the MLR-estimated value.

        Physical Similarity with IDW and regression-based augmentation (PS_IDW_RA)
            The ungauged hydrograph is an average of the `n` most similar catchments' hydrographs, but
            the average is weighted using inverse distance weighting. Furthermore, the method uses the CANOPEX/USGS
            dataset to estimate model parameters using Multiple Linear Regression. Parameters whose regression r-squared
            is higher than 0.5 are replaced by the MLR-estimated value.

    """
    identifier = "regionalisation"
    title = "Simulate streamflow at an ungauged site based on surrounding or similar gauged catchments."
    abstract = "Compute the hydrograph for an ungauged catchment using a regionalization method."
    version = '0.1'

    method = LiteralInput('method', 'Regionalisation method',
                          abstract="Regionalisation method to use, one of MLR, SP, PS, SP_IDW, "
                                   "PS_IDW, SP_IDW_RA, PS_IDW_RA.",
                          data_type='string',
                          allowed_values=(
                              'MLR', 'SP', 'PS', 'SP_IDW', 'PS_IDW', 'SP_IDW_RA', 'PS_IDW_RA'),
                          default='SP_IDW',
                          min_occurs=0)

    ndonors = LiteralInput('ndonors', 'Number of gauged catchments to use for the regionalizaion.',
                           abstract="Number of close or similar catchments to use to generate the representative "
                                    "hydrograph at the ungauged site.",
                           data_type='integer',
                           default=5,
                           min_occurs=0)

    min_NSE = LiteralInput('min_NSE', 'NSE Score (unitless)',
                           abstract="Minimum calibration NSE value required to be considered in the regionalization.",
                           data_type='float',
                           default=0.6,
                           min_occurs=0)

    properties = ComplexInput('properties', 'Regionalization properties',
                              abstract="json string storing dictionary of properties.",
                              min_occurs=1,
                              max_occurs=1,
                              supported_formats=[FORMATS.JSON, ])

    inputs = [wio.ts, wio.start_date, wio.end_date, wio.latitude, wio.longitude,
              wio.model_name, ndonors, min_NSE, method, properties, wio.area, wio.elevation]

    outputs = [wio.hydrograph, wio.ensemble]

    def _handler(self, request, response):
        response.update_status('PyWPS process {} started.'.format(self.identifier), 0)

        ts = [e.file for e in request.inputs.pop('ts')]
        model_name = request.inputs.pop('model_name')[0].data
        method = request.inputs.pop('method')[0].data
        ndonors = request.inputs.pop('ndonors')[0].data
        latitude = request.inputs.pop('latitude')[0].data
        longitude = request.inputs.pop('longitude')[0].data
        min_NSE = request.inputs.pop('min_NSE')[0].data
        properties = request.inputs.pop('properties')[0].data
        properties = json.loads(properties)

        kwds = {}
        for key, val in request.inputs.items():
            kwds[key] = request.inputs[key][0].data

        nash, params = read_gauged_params(model_name)
        props = read_gauged_properties(properties)

        ungauged_props = {key: properties[key] for key in properties}
        # kwds.update(properties) # This fails as properties are not part of the Raven keywords (i.e. "forest")

        qsim, ensemble = regionalize(method, model_name, nash, params,
                                     props, ungauged_props,
                                     size=ndonors,
                                     min_NSE=min_NSE,
                                     ts=ts,
                                     **kwds)

        # Write output
        nc_qsim = Path(self.workdir) / 'qsim.nc'
        qsim.to_netcdf(nc_qsim)
        response.outputs['hydrograph'].file = str(nc_qsim)

        # TODO: Commplete attributes
        nc_ensemble = Path(self.workdir) / 'ensemble.nc'
        ensemble.to_netcdf(nc_ensemble)
        response.outputs['ensemble'].file = str(nc_ensemble)

        return response
