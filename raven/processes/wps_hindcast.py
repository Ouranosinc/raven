import logging

from pathlib import Path
import fiona
from pywps import ComplexInput, LiteralInput, Process, FORMATS
import tempfile
import pandas as pd

from raven.utilities import forecasting
from raven.utils import archive_sniffer, single_file_check

from . import wpsio as wio

LOGGER = logging.getLogger("PYWPS")

class HindcastingProcess(Process):
    def __init__(self):
        """
        Notes
        -----
        This WPS service provides hindcasting capabilities for a given date in
        the ECCC forecasting archive as made available by the CaSPAr (Canadian
        Surface Prediction Archive) as can be found here: https://caspar-data.ca/
        with more information here:
        
            Mai, J., Kornelsen, K.C., Tolson, B.A., Fortin, V., Gasset, N., Bouhemhem, D.,
        Schäfer, D., Leahy, M., Anctil, F. and Coulibaly, P., 2020. The Canadian 
        Surface Prediction Archive (CaSPAr): A Platform to Enhance Environmental 
        Modeling in Canada and Globally. Bulletin of the American Meteorological 
        Society, 101(3), pp.E341-E356.
        """
            
        climate_model = LiteralInput('climate_model', 'Which ECCC forecast model should be used',
                              abstract='Which of the four climate models should be used: GEPS, GDPS, REPS, RDPS.',
                              data_type='string',
                              allowed_values=('GEPS', 'REPS', 'GDPS', 'RDPS'),
                              default='GEPS',
                              min_occurs=1)
        
        region_vector = ComplexInput('region_vector', 'Vector Shape of the desired location',
                             abstract='An ESRI Shapefile, GML, JSON, GeoJSON, or single layer GeoPackage.'
                                      ' The ESRI Shapefile must be zipped and contain the .shp, .shx, and .dbf.',
                             min_occurs=1, max_occurs=1,
                             supported_formats=[FORMATS.GEOJSON, FORMATS.GML, FORMATS.JSON, FORMATS.SHP, FORMATS.ZIP])
        
        params = LiteralInput("params","Comma separated list of model parameters",
                              abstract="Parameters to run the model",
                              data_type="string", min_occurs=1, max_occurs=1)

        hdate =  LiteralInput('hdate', 'Hindcast start date (AAAA-MM-DD)',
                          abstract='Start date of the hindcast (AAAA-MM-DD). '
                                   'Defaults to the start of the forcing file. ',
                          data_type='dateTime', min_occurs=1, max_occurs=1)
                             
        inputs = [params, wio.latitude, wio.longitude, wio.name, wio.area, wio.elevation, 
                  wio.model_name, climate_model, region_vector, wio.rain_snow_fraction, wio.nc_spec, wio.rvc, hdate]
    
        outputs = [wio.forecast]
    
        super(HindcastingProcess, self).__init__(
                self._handler,
                identifier = "hindcasting",
                title = "Perform hindcast forecast for a given forecast date in CaSPAr.",
                abstract = "Perform a deterministic or probabilistic raven forecast using a given forecast date in CaSPAr.",
                version = '0.1',
                inputs=inputs,
                outputs=outputs,
                keywords=["hindcasting", "CaSPAr", "GEPS", "REPS", "GDPS", "RDPS", "ensemble forecasts"],
                status_supported=True,
                store_supported=True)
        
    # Start WPS
    def _handler(self, request, response):
        
       
        # Get what is needed integrally
        model_name = request.inputs.pop('model_name')[0].data
        climate_model = request.inputs.pop('climate_model')[0].data
        region_vector = request.inputs.pop('region_vector')[0].file
        rvc=request.inputs.pop('rvc')[0].file
        hdate=request.inputs.pop('hdate')[0].data
        
        extensions = ['.gml', '.shp', '.gpkg', '.geojson', '.json']
        vector_file = single_file_check(archive_sniffer(region_vector, working_dir=self.workdir, extensions=extensions))
        
        
        # Build the other kewds needed to run the model
        kwds = {}
        for key, val in request.inputs.items():
            kwds[key] = request.inputs[key][0].data
        response.update_status('Inputs are read', 1)
        
        
        # Prepare the forecast data from the latest ECCC forecast.
        hcst = forecasting.get_hindcast_day(
               fiona.open(vector_file), hdate, climate_model=climate_model)
        
        response.update_status('Hindcast data is collected', 2)

        # write the forecast data to file
        tmp_dir = tempfile.TemporaryDirectory()
        hcst.to_netcdf(f"{tmp_dir.name}/hcstfile.nc")
        
        # Add the forecast file to the ts needed by Raven and run the model
        kwds['ts']=f"{tmp_dir.name}/hcstfile.nc"
        kwds['nc_index']=range(hcst.dims.get("member"))
        kwds['start_date']=pd.to_datetime(hcst.time[0].values)
      
        qsim=forecasting.perform_forecasting_step(rvc, model_name, **kwds)
        
        response.update_status('Computed hydrological forecast', 99)

        # Write output
        nc_qsim = Path(self.workdir) / 'qsim.nc'
        qsim.to_netcdf(nc_qsim)
        response.outputs['forecast'].file = str(nc_qsim)

        return response