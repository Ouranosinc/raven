import logging
import json
import tempfile

from pywps import LiteralInput, ComplexInput
from pywps import ComplexOutput
from pywps import Process, FORMATS
from pywps.app.Common import Metadata

from raven.utils import archive_sniffer, crs_sniffer, single_file_check, generic_vector_reproject
from raven.utilities import gis
import pdb

LOGGER = logging.getLogger("PYWPS")



class GetGeometDataProcess(Process):
    """Extract GeoMet forecast data for doing hydrological forecasting for a watershed whose shape is passed as a parameter"""

    def __init__(self):
        inputs = [
            ComplexInput('shape', 'Vector Shape',
                         abstract='An ESRI Shapefile, GML, JSON, GeoJSON, or single layer GeoPackage.'
                                  ' The ESRI Shapefile must be zipped and contain the .shp, .shx, and .dbf.'
                                  ' The shape and raster should have a matching CRS.',
                         min_occurs=1, max_occurs=1,
                         supported_formats=[FORMATS.GEOJSON, FORMATS.GML, FORMATS.JSON, FORMATS.SHP]),
            LiteralInput('forecast_date', 'which forecast date should be extracted?',
                         data_type='dateTime',
                         abstract='provide the date of the forecast issue',
                         min_occurs=1, max_occurs=1),
            LiteralInput('use_rdps', 'Should the forecast include the RDPS (regional deterministic) model?',
                         data_type='boolean', default=True,
                         abstract='Should the RDPS model be extracted?',
                         min_occurs=0, max_occurs=1),
            LiteralInput('use_gdps', 'Should the forecast include the RDPS (global deterministic) model?',
                         data_type='boolean', default=True,
                         abstract='Should the GDPS model be extracted?',
                         min_occurs=0, max_occurs=1),
            LiteralInput('use_reps', 'Should the forecast include the REPS (regional ensemble) model?',
                         data_type='boolean', default=True,
                         abstract='Should the REPS model be extracted?',
                         min_occurs=0, max_occurs=1),
            LiteralInput('use_geps', 'Should the forecast include the GEPS (global ensemble) model?',
                         data_type='boolean', default=True,
                         abstract='Should the GEPS model be extracted?',
                         min_occurs=0, max_occurs=1),
            LiteralInput('combine_reg_glob', 'agregate regional and global?',
                         data_type='boolean', default=True,
                         abstract='Should the forecasts be agregated with REG + global after regional model forecast ends?',
                         min_occurs=0, max_occurs=1),
        ]

        outputs = [
            ComplexOutput('rdps', 'RDPS forecast netcdf file averaged at the catchment scale',
                          abstract='RDPS forecast netcdf file averaged at the catchment scale',
                          supported_formats=[FORMATS.NETCDF]),
            ComplexOutput('gdps', 'GDPS forecast netcdf file averaged at the catchment scale',
                          abstract='GDPS forecast netcdf file averaged at the catchment scale',
                          supported_formats=[FORMATS.NETCDF]),
            ComplexOutput('reps', 'REPS forecast netcdf file averaged at the catchment scale',
                          abstract='REPS forecast netcdf file averaged at the catchment scale',
                          supported_formats=[FORMATS.NETCDF]),
            ComplexOutput('geps', 'GEPS forecast netcdf file averaged at the catchment scale',
                          abstract='GEPS forecast netcdf file averaged at the catchment scale',
                          supported_formats=[FORMATS.NETCDF]),
            ComplexOutput('rdps_gdps', 'combined RDPS and GDPS forecast netcdf file averaged at the catchment scale',
                          abstract='combined RDPS and GDPS forecast netcdf file averaged at the catchment scale',
                          supported_formats=[FORMATS.NETCDF]),
            ComplexOutput('reps_geps', 'combined REPS and GEPS forecast netcdf file averaged at the catchment scale',
                          abstract='combined REPS and GEPS forecast netcdf file averaged at the catchment scale',
                          supported_formats=[FORMATS.NETCDF]),
            
        ]

        super(GetGeometDataProcess, self).__init__(
            self._handler,
            identifier="geomet-forecasts",
            title="Get Geomet forecast data",
            version="1.0",
            abstract="Get geomet forecast data for a given day to drive the RAVEN models",
            metadata=[],
            inputs=inputs,
            outputs=outputs,
            status_supported=True,
            store_supported=True)

    def _handler(self, request, response):

        # Variables that we need for RAVEN:
        GeoMetVariables=['ETA_PR','ETA_TT','ETA_RN','ETA_SN'] # In order: Precip, air temp, rain and snow
        shape_url = request.inputs['shape'][0].file
        use_rdps = request.inputs['use_rdps'][0].data
        use_gdps = request.inputs['use_gdps'][0].data
        use_reps = request.inputs['use_reps'][0].data
        use_geps = request.inputs['use_geps'][0].data
        combine_reg_glob = request.inputs['combine_reg_glob'][0].data
        forecast_date = request.inputs['forecast_date'][0].data
        vectors = ['.gml', '.shp', '.gpkg', '.geojson', '.json']
        vector_file = single_file_check(archive_sniffer(shape_url, working_dir=self.workdir, extensions=vectors))
        
        bbox = gis.get_bbox(vector_file)

        models=[]
        if use_rdps:
            models.append("RDPS")
        if use_gdps:
            models.append("GDPS")
        if use_reps:
            models.append("REPS")
        if use_geps:
            models.append("GEPS")
        pdb.set_trace()    
                
        for MODEL in models:
            for VARIABLE in GeoMetVariables:
                # Careful: Time is finicky, need to make sure the date is ok: Ex: on 2020-04-26, 21:20, the first available forecast was for 2020-04-27 T18:00:00Z. We can find the latest forecast run using a WMS GetCapabilities request, to do later.
                query='https://geo.weather.gc.ca/geomet?SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage&COVERAGEID=' + MODEL + '.' + VARIABLE + '&SUBSETTINGCRS=EPSG:4326&SUBSET=x(' + str(bbox[0]-0.5) + ',' + str(bbox[2]+0.5) + ')&SUBSET=y(' + str(bbox[1]-0.5) + ',' + str(bbox[3]+0.5) + ')&FORMAT=image/netcdf&TIME=' + forecast_date.strftime("%Y-%m-%d") + 'T06:00:00Z'
        
            
        response.outputs['rdps'].data = 1


        return response
