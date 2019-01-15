# -*- coding: utf-8 -*-


import datetime as dt


from pywps import Service
from pywps.tests import assert_response_success
from raven.processes import RegionalisationProcess
from tests.common import client_for, TESTDATA, get_output, urlretrieve
from netCDF4 import Dataset
import pdb


# @pytest.mark.skip
class TestRegionalisation:
    
  
    
    def testRegionalisationHMETS_1(self):
         
        number_donors=1
        

        client = client_for(Service(processes=[RegionalisationProcess(),]))
             
        datainputs = "ts=files@xlink:href=file://{ts};" \
                     "start_date={start_date};" \
                     "end_date={end_date};" \
                     "name={name};" \
                     "area={area};" \
                     "latitude={latitude};" \
                     "longitude={longitude};" \
                     "elevation={elevation};" \
                     "model_name={model_name};" \
                     "min_NSE={min_NSE};" \
                     "number_donors={number_donors};" \
                     "regionalisationMethod={regionalisationMethod};"\
                .format(ts=TESTDATA['raven-hmets-nc-ts'],
                    start_date=dt.datetime(2000, 1, 1),
                    end_date=dt.datetime(2002, 1, 1),
                    name='Salmon',
                    run_name='test',
                    area='4250.6',
                    elevation='843.0',
                    latitude=40.4848,
                    longitude=-103.3659,
                    model_name='HMETS',
                    min_NSE=0.6,
                    number_donors=number_donors,
                    regionalisationMethod='PS_IDW_RA',
                    )
             
         
        resp = client.get(service='WPS', request='execute', version='1.0.0', identifier='RegionalisationProcess', datainputs=datainputs)
         
         
        assert_response_success(resp)
         
         
        out = get_output(resp.xml)
        assert 'hydrograph' in out
        tmp_file, _ = urlretrieve(out['hydrograph'])
        nc_file=Dataset(tmp_file)
        assert(nc_file.variables['Qsim'].shape[1]==number_donors)
             
   
    def testRegionalisationHMETS_2(self):
         
        number_donors=5
        
    
            
        client = client_for(Service(processes=[RegionalisationProcess(),]))
    
        datainputs = "ts=files@xlink:href=file://{ts};" \
                     "start_date={start_date};" \
                     "end_date={end_date};" \
                     "name={name};" \
                     "area={area};" \
                     "latitude={latitude};" \
                     "longitude={longitude};" \
                     "elevation={elevation};" \
                     "model_name={model_name};" \
                     "min_NSE={min_NSE};" \
                     "number_donors={number_donors};" \
                     "regionalisationMethod={regionalisationMethod};"\
                .format(ts=TESTDATA['raven-hmets-nc-ts'],
                    start_date=dt.datetime(2000, 1, 1),
                    end_date=dt.datetime(2002, 1, 1),
                    name='Salmon',
                    run_name='test',
                    area='4250.6',
                    elevation='843.0',
                    latitude=54.4848,
                    longitude=-123.3659,
                    model_name='HMETS',
                    min_NSE=0.6,
                    number_donors=number_donors,
                    regionalisationMethod='PS_IDW',
                    )

    
         
        resp = client.get(service='WPS', request='execute', version='1.0.0', identifier='RegionalisationProcess', datainputs=datainputs)
         
         
        assert_response_success(resp)
         
         
        out = get_output(resp.xml)
        assert 'hydrograph' in out
        tmp_file, _ = urlretrieve(out['hydrograph'])
        nc_file=Dataset(tmp_file)
        assert(nc_file.variables['Qsim'].shape[1]==number_donors)




    def testRegionalisationGR4J_1(self):
         
            number_donors=1

            client = client_for(Service(processes=[RegionalisationProcess(),]))
    
            datainputs = "ts=files@xlink:href=file://{ts};" \
                     "start_date={start_date};" \
                     "end_date={end_date};" \
                     "name={name};" \
                     "area={area};" \
                     "latitude={latitude};" \
                     "longitude={longitude};" \
                     "elevation={elevation};" \
                     "model_name={model_name};" \
                     "min_NSE={min_NSE};" \
                     "number_donors={number_donors};" \
                     "regionalisationMethod={regionalisationMethod};"\
                    .format(ts=TESTDATA['raven-gr4j-cemaneige-nc-ts'],
                    start_date=dt.datetime(2000, 1, 1),
                    end_date=dt.datetime(2002, 1, 1),
                    name='Salmon',
                    run_name='test',
                    area='4250.6',
                    elevation='843.0',
                    latitude=40.4848,
                    longitude=-105.3659,
                    model_name='GR4JCN',
                    min_NSE=0.8,
                    number_donors=number_donors,
                    regionalisationMethod='SP_IDW',
                    )

         
            resp = client.get(service='WPS', request='execute', version='1.0.0', identifier='RegionalisationProcess', datainputs=datainputs)
         
         
            assert_response_success(resp)
         
         
            out = get_output(resp.xml)
            assert 'hydrograph' in out
            tmp_file, _ = urlretrieve(out['hydrograph'])
            nc_file=Dataset(tmp_file)
            assert(nc_file.variables['Qsim'].shape[1]==number_donors)
             
        

    def testRegionalisationGR4J_2(self):
         
             number_donors=5
        
             client = client_for(Service(processes=[RegionalisationProcess(),]))
    
             datainputs = "ts=files@xlink:href=file://{ts};" \
                     "start_date={start_date};" \
                     "end_date={end_date};" \
                     "name={name};" \
                     "area={area};" \
                     "latitude={latitude};" \
                     "longitude={longitude};" \
                     "elevation={elevation};" \
                     "model_name={model_name};" \
                     "min_NSE={min_NSE};" \
                     "number_donors={number_donors};" \
                     "regionalisationMethod={regionalisationMethod};"\
                    .format(ts=TESTDATA['raven-gr4j-cemaneige-nc-ts'],
                    start_date=dt.datetime(2000, 1, 1),
                    end_date=dt.datetime(2002, 1, 1),
                    name='Salmon',
                    run_name='test',
                    area='4250.6',
                    elevation='843.0',
                    latitude=40.4848,
                    longitude=-105.3659,
                    model_name='GR4JCN',
                    min_NSE=0.6,
                    number_donors=number_donors,
                    regionalisationMethod='SP_IDW_RA',
                    )

    
         
             resp = client.get(service='WPS', request='execute', version='1.0.0', identifier='RegionalisationProcess', datainputs=datainputs)
         
         
             assert_response_success(resp)
         
         
             out = get_output(resp.xml)
             assert 'hydrograph' in out
             tmp_file, _ = urlretrieve(out['hydrograph'])
             nc_file=Dataset(tmp_file)
             assert(nc_file.variables['Qsim'].shape[1]==number_donors)
                          
             
    def testRegionalisationMOHYSE_1(self):
         
             number_donors=1
        
        

             client = client_for(Service(processes=[RegionalisationProcess(),]))
    
             datainputs = "ts=files@xlink:href=file://{ts};" \
                     "start_date={start_date};" \
                     "end_date={end_date};" \
                     "name={name};" \
                     "area={area};" \
                     "latitude={latitude};" \
                     "longitude={longitude};" \
                     "elevation={elevation};" \
                     "model_name={model_name};" \
                     "min_NSE={min_NSE};" \
                     "number_donors={number_donors};" \
                     "regionalisationMethod={regionalisationMethod};"\
                    .format(ts=TESTDATA['raven-mohyse-nc-ts'],
                    start_date=dt.datetime(2000, 1, 1),
                    end_date=dt.datetime(2002, 1, 1),
                    name='Salmon',
                    run_name='test',
                    area='4250.6',
                    elevation='843.0',
                    latitude=40.4848,
                    longitude=-105.3659,
                    model_name='MOHYSE',
                    min_NSE=0.7,
                    number_donors=number_donors,
                    regionalisationMethod='SP',
                    )

         
             resp = client.get(service='WPS', request='execute', version='1.0.0', identifier='RegionalisationProcess', datainputs=datainputs)
         
         
             assert_response_success(resp)
         
         
             out = get_output(resp.xml)
             assert 'hydrograph' in out
             tmp_file, _ = urlretrieve(out['hydrograph'])
             nc_file=Dataset(tmp_file)
             assert(nc_file.variables['Qsim'].shape[1]==number_donors)
             
        

    def testRegionalisationMOHYSE2_2(self):
         
             number_donors=5
         
             client = client_for(Service(processes=[RegionalisationProcess(),]))
    
             datainputs = "ts=files@xlink:href=file://{ts};" \
                     "start_date={start_date};" \
                     "end_date={end_date};" \
                     "name={name};" \
                     "area={area};" \
                     "latitude={latitude};" \
                     "longitude={longitude};" \
                     "elevation={elevation};" \
                     "model_name={model_name};" \
                     "min_NSE={min_NSE};" \
                     "number_donors={number_donors};" \
                     "regionalisationMethod={regionalisationMethod};"\
                    .format(ts=TESTDATA['raven-mohyse-nc-ts'],
                    start_date=dt.datetime(2000, 1, 1),
                    end_date=dt.datetime(2002, 1, 1),
                    name='Salmon',
                    run_name='test',
                    area='4250.6',
                    elevation='843.0',
                    latitude=40.4848,
                    longitude=-105.3659,
                    model_name='MOHYSE',
                    min_NSE=0.6,
                    number_donors=number_donors,
                    regionalisationMethod='PS',
                    )

    
         
             resp = client.get(service='WPS', request='execute', version='1.0.0', identifier='RegionalisationProcess', datainputs=datainputs)
         
         
             assert_response_success(resp)
         
         
             out = get_output(resp.xml)
             assert 'hydrograph' in out
             tmp_file, _ = urlretrieve(out['hydrograph'])
             nc_file=Dataset(tmp_file)
             assert(nc_file.variables['Qsim'].shape[1]==number_donors)
    