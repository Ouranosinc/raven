
"""
entrypoint_regionalization.py

INPUTS:
    hydrological_model -- The hydrological model to use, limited as of now to:
                          'HMETS', 'MOHYSE', 'GR4J'.
    inputs_file -- Ungauged catchment inputs file containing the following fields:
        rain (extracted from boundaries beforehand)
        snow (extracted from boundaries beforehand)
        tasmax (extracted from boundaries beforehand)
        tasmin (extracted from boundaries beforehand)
        surface_area (extracted from GIS beforehand)
        day of year (extracted from dataset used for meteo)
        dates (extracted from dataset used for meteo)
        latitude (centroid, extracted from GIS beforehand)
        longitude (centroid, extracted from GIS beforehand)
        land cover fractions (determined by GIS tools beforehand)
        avg slope (determined by GIS beforehand)
        avg elevation (determined by GIS beforehand)

    regionalization_method -- An integer with the value-argument pairs:
        1- Multiple Linear Regression (MLR)
        2- Spatial Proximity (SP)
        3- Spatial Proximity with Inverse Distance Weighting (IDW) of
           donor basins (SP_IDW)
        4- Spatial Proximity with IDW and regression-based augmentation (RA)
            (SP_IDW_RA)
        5- Physical Similarity (PS)
        6- Physical Similarity with IDW of donor basins (PS_IDW)
        7- Physical Similarity with IDW and RA (PS_IDW_RA)

    number_donors -- A positive integer from [1:10] that allows averaging the
        hydrographs of number_donors donors. For example, a value of 5
        indicates that we take the average hydrograph of the 5 closest or most
        similar donors. No impact when used for MLR. IDW computes the average 
        using inverse weighting.

    min_NSE -- The minimum calibration NSE required for a donor catchment to be
        considered. If the closest or most similar catchment has an NSE below
        the min_NSE, then it is skipped and the next closest/similar catchment
        is selected.

OUTPUT:
    This function returns two arrays. Qsim = [n_days x number_donors] 
    in which the columns are the individual hydrographs used to produce the 
    averaged hydrograph, in order of closest to farthest (or most to least similar)
    The second is "best_estimate", which is the average (or IDW average) of the 
    Qsim hydrographs.
"""


from . import wpsio as wio
import logging
from pywps import Process, LiteralInput
from pathlib import Path
from . regionalization_tools import regionalization
from .wps_raven import RavenProcess
LOGGER = logging.getLogger("PYWPS")


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
            ...

        Physical similarity with inverse distance weighting (PS_IDW)
            ...

        Spatial proximity with IDW and regression-based augmentation (SP_IDW_RA)
            ...

        Physical Similarity with IDW and regression-based augmentation (PS_IDW_RA)
            ...

    """
    identifier = "regionalisation"
    title = "Simulate streamflow at an ungauged site based on surrounding or similar gauged catchments."
    abstract = "Compute the hydrograph for an ungauged catchment using a regionalization method."
    version = '0.1'

    method = LiteralInput('method', 'Regionalisation method',
                                          abstract="Regionalisation method to use, one of MLR, SP, PS, SP_IDW, "
                                                   "PS_IDW, SP_IDW_RA,"
                                                   "PS_IDW_RA.",
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

    inputs = [wio.ts, wio.start_date, wio.end_date, wio.area, wio.elevation, wio.latitude, wio.longitude,
              wio.model_name, ndonors, min_NSE, method]

    outputs = [wio.hydrograph, wio.ensemble]

    def _handler(self, request, response):
        response.update_status('PyWPS process {} started.'.format(self.identifier), 0)

        keys = request.inputs.keys()

        ts = [e.file for e in keys.pop('ts')]
        model_name = keys.pop('model_name')[0].data
        method = keys.pop('method')[0].data
        ndonors = keys.pop('ndonors')[0].data
        latitude = keys.pop('latitude')[0].data
        longitude = keys.pop('longitude')[0].data
        min_NSE = keys.pop('min_NSE')[0].data

        kwds = {}
        for key in keys:
            kwds[key] = request.inputs[key][0].data

        qsim, ensemble = regionalization(method, model_name,
                                         latitude=latitude, longitude=longitude,
                                         size=ndonors,
                                         min_NSE=min_NSE,
                                         **kwds)

        # Write output
        nc_qsim = Path(self.workdir) / 'qsim.nc'
        qsim.to_dataset(nc_qsim)
        response.ouputs['hydrograph'].file = nc_qsim

        nc_ensemble = Path(self.workdir) / 'ensemble.nc'
        qsim.to_dataset(nc_ensemble)
        response.ouputs['ensemble'].file = nc_ensemble

        return response
