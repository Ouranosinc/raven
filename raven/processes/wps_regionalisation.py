import logging
from pathlib import Path

import xarray as xr
from pywps import Process, LiteralInput

from raven.utilities import regionalize, read_gauged_properties, read_gauged_params
from . import wpsio as wio
from .wps_raven import RavenProcess

LOGGER = logging.getLogger("PYWPS")


# TODO: latitude and longitude have a different meaning here if we're using them to get the catchment properties.
#  Normally for other WPS Raven processes, they refer to the centroid, here they'd refer to the outlet, correct ?
#  ANSWER: No, we're still talking about the centroid! basically the closer the center of mass of the catchment, the
#          more likely the catchments will be physically and hydrologically similar, according to the philosophy
#          behind the method.
# But I mean in this process, aren't we passing lat, lon to extract properties for a new watershed? And then we'll
# extract the centroid lat and lon for the analysis.


class RegionalisationProcess(RavenProcess):
    """
    TODO: Include a description of each method.
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

    inputs = [wio.ts, wio.start_date, wio.end_date, wio.latitude, wio.longitude,
              wio.model_name, ndonors, min_NSE, method]

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

        kwds = {}
        for key, val in request.inputs.items():
            kwds[key] = request.inputs[key][0].data

        nash, params = read_gauged_params(model_name)
        variables = ['longitude', 'latitude']
        props = read_gauged_properties()[variables]

        # TODO: Replace by function determining catchment properties from DEM and land use file and Hydrosheds data.
        def get_catchment_properties(lat, lon):
            return {'longitude': .7, 'latitude': .7, 'area': '4250.6', 'elevation': '843.0'}

        catchment_props = get_catchment_properties(latitude, longitude)
        properties = ['longitude', 'latitude']
        ungauged_props = {key: catchment_props[key] for key in properties}
        kwds.update(catchment_props)

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
