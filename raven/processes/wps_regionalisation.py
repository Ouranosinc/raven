
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
from pywps import Process
import numpy as np
from netCDF4 import Dataset

LOGGER = logging.getLogger("PYWPS")

from . regionalization_tools import (
        get_array_of_gauged_properties,
        filter_array_by_NSE,
        rank_donor_catchments,
        get_ungauged_properties,
        IDW
)

from . regionalization_methods import (
        MLR_parameters_regionalization,
        distance_only_regionalization,
        distance_MLR_regionalization
)

import pdb

class RegionalisationProcess(Process):
    identifier = 'RegionalisationProcess'
    abstract = 'Regionalisation methods to predict streamflow at ungauged locations'
    title = "Process to predict streamflow at ungauged sites based on surrounding or similar gauged catchments"
    version = '0.1'

    inputs = [wio.ts, wio.start_date, wio.end_date, wio.area, wio.elevation, wio.latitude, wio.longitude, wio.model_name, wio.number_donors, wio.min_NSE, wio.regionalisationMethod]

    outputs = [wio.hydrograph]

    

    def __init__(self):
       


        super(RegionalisationProcess, self).__init__(
            self._handler,
            identifier="RegionalisationProcess",
            title="Methods to predict streamflow at ungauged sites",
            version="0.1",
            abstract="Use this method to generate streamflow at ungauged sites using a rich database of catchments over North America.",
            metadata=[],
            inputs=self.inputs,
            outputs=self.outputs,
            status_supported=True,
            store_supported=True)
       
        

    def _handler(self, request, response):

       
        
        response.update_status('PyWPS process {} started.'.format(self.identifier), 0)
        
      
        # Model name
        if 'model_name' in request.inputs:
            model_name=request.inputs['model_name'][0].data
        else:
            raise ValueError('Model Name is not optional')
       
        #Regionalisation Method
         
        if 'regionalisationMethod' in request.inputs:
            regionalisationMethod=request.inputs['regionalisationMethod'][0].data
        else:
            regionalisationMethod=2
            print('regionalisationMethod not provided, set to 2 by default')
   
        print("REGIONALIZATION METHOD:")
        print(regionalisationMethod)
        print(" IN PROGRESS")
             
        # Number Donors
        if 'number_donors' in request.inputs:
            number_donors=request.inputs['number_donors'][0].data
        else:
            number_donors=5
            print('Number of donors not provided, set to 5 by default')

        # Minimum NSE filter
        if 'min_NSE' in request.inputs:
            min_NSE=request.inputs['min_NSE'][0].data
        else:
            min_NSE=0.6
            print('minimum NSE not provided, set to 0.6 by default')

        # Model inputs
        ts=[f.file for f in request.inputs['ts']]
        
        
        ncfilename=request.inputs['ts'][0].file
        start_date=request.inputs['start_date'][0].data
        end_date=request.inputs['end_date'][0].data
        latitude=request.inputs['latitude'][0].data
        longitude=request.inputs['longitude'][0].data
        elevation=request.inputs['elevation'][0].data
        area=request.inputs['area'][0].data
        
        # Find the properties we wish to use for multiple linear regression (MLR)
        # and similarity.
        """
        properties_to_use = ['centroid_latitude',
                         'centroid_longitude',
                         'surface_area',
                         'avg_elevation',
                         'avg_slope',
                         'land_forest',
                         'land_grass',
                         'land_impervious',
                         'land_urban']
        """
        properties_to_use=['centroid_latitude','centroid_longitude']
        # Get the ungauged catchment properties from the inputs_file for the
        # regionalization scheme.
        ungauged_properties = get_ungauged_properties(ts,properties_to_use,latitude,longitude)
        

        # Load gauged catchments properties and hydrological model parameters.
        gauged_properties, model_parameters = get_array_of_gauged_properties(
                                            model_name)

        # Filter arrays according to the minimum value allowed.
        gauged_properties, model_parameters = filter_array_by_NSE(
                                            gauged_properties,
                                            model_parameters,
                                            min_NSE)

        # Check to see if we have enough data, otherwise raise error
        if (gauged_properties.shape[0] < number_donors and
                regionalisationMethod >= 2):
            raise ValueError("hydrological_model and minimum NSE threshold \
                             combination is too strict for the number of donor \
                             basins. Please reduce the number of donor basins OR \
                             reduce the minimum NSE threshold")

        # Rank the matrix according to the similarity or distance.
        if (regionalisationMethod >= 5):  # physical similarity
            regionalisation_tag = 'similarity'
        else:
            # Filter by distance if MLR (option 1) is selected, it is faster
            # than by similarity.
            regionalisation_tag = 'distance'

            # Create a dataframe that is identical to gauged_properties arrays.
        gauged_properties, model_parameters = rank_donor_catchments(
                                              ungauged_properties,
                                              model_parameters,
                                              gauged_properties,
                                              properties_to_use,
                                              regionalisation_tag
                                              )

        if regionalisationMethod == 1:
            # Multiple linear regression (MLR)
            Qsim = MLR_parameters_regionalization(ncfilename,start_date, end_date, latitude, longitude, elevation, area,
                                                  model_name, ts,
                                                  ungauged_properties,
                                                  model_parameters,
                                                  gauged_properties,
                                                  properties_to_use)
            best_estimate = Qsim

        elif (regionalisationMethod == 2) or (regionalisationMethod == 5):
            # Simple averaging
            Qsim = distance_only_regionalization(ncfilename,start_date, end_date, latitude, longitude, elevation, area,
                                                 model_name,
                                                 ts,
                                                 number_donors,
                                                 model_parameters)
            best_estimate = np.average(Qsim, 1)

        elif (regionalisationMethod == 3) or (regionalisationMethod == 6):
            # Inverse distance weighting (IDW)
            Qsim = distance_only_regionalization(ncfilename,start_date, end_date, latitude, longitude, elevation, area,
                                                 model_name,
                                                 ts,
                                                 number_donors,
                                                 model_parameters)
         
            
            best_estimate = IDW(Qsim, number_donors, gauged_properties['distance'])

        elif (regionalisationMethod == 4) or (regionalisationMethod == 7):
            # Inverse distance weighting (IDW)
            Qsim = distance_MLR_regionalization(ncfilename,start_date, end_date, latitude, longitude, elevation, area,
                                                model_name, ts,
                                                number_donors,
                                                ungauged_properties,
                                                model_parameters,
                                                gauged_properties,
                                                properties_to_use)
            
            best_estimate = IDW(Qsim, number_donors, gauged_properties['distance'])

        else:
            print("Unknown regionalization method")


        # Write utput
        indexfile = 'RegionalisationHydrograph.nc'
       # time = pd.date_range(start='2000-07-01', end='2002-07-01', freq='D')
     
        with Dataset(indexfile, 'w', format='NETCDF4') as nc_file:
            nc_file.createDimension('time', Qsim.shape[0])
            nc_file.createDimension('donor_number', Qsim.shape[1])
          #  nc_file.createDimension('AveragedEstimate',1)
            vars_list = nc_file.createVariable('Qsim', 'float32', ('time', 'donor_number'))
            vars_list[:] = Qsim
            vars_list2=nc_file.createVariable('BestEstimate','float32',('time'))
            vars_list2[:]=best_estimate
           # vars_list3=nc_file.createVariable('Time','datetime64',('time'))
           # vars_list3[:]=time
        #response.outputs['hy'] = Qsim
        
        
        response.outputs['hydrograph'].file=indexfile
        return response
    