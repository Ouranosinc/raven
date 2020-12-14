import pytest
from pywps import Service
from pywps.tests import assert_response_success


from raven.processes import GraphFcstUncertaintyProcess
from .common import client_for, TESTDATA, CFG_FILE, get_output


class TestForecastGraphProcess:

    def test_forecast_graph(self):
        client = client_for(Service(processes=[GraphFcstUncertaintyProcess(), ], cfgfiles=CFG_FILE))
        
        datainputs = "fcst=files@xlink:href=file://{fcst};" \
                     "fcst_var={fcst_var};".format(fcst=TESTDATA['floodrisk_ens'], fcst_var='fcst')
  
        resp = client.get(service='WPS', request='Execute', version='1.0.0', identifier='graph_forecast_uncertainty',datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)
        assert out['graph_forecasts'].endswith('.png')
        
    def test_hindcast_graph(self):
        client = client_for(Service(processes=[GraphFcstUncertaintyProcess(), ], cfgfiles=CFG_FILE))
        
        
        datainputs = "fcst=files@xlink:href=file://{fcst};" \
                     "qobs=files@xlink:href=file://{qobs};" \
                     "fcst_var={fcst_var};" \
                     "qobs_var={qobs_var};".format(fcst=TESTDATA['floodrisk_ens'], fcst_var='fcst', qobs=TESTDATA['XSS_obs'], qobs_var='obs')
  
        resp = client.get(service='WPS', request='Execute', version='1.0.0', identifier='graph_forecast_uncertainty',datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)
        assert out['graph_forecasts'].endswith('.png')
        
     