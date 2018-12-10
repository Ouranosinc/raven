import os
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
from pywps.tests import WpsClient, WpsTestResponse
from pywps import get_ElementMakerForVersion
from pywps.app.basic import get_xpath_ns
import six
if six.PY2:
    from urllib import urlretrieve
else:
    from urllib.request import urlretrieve


VERSION = "1.0.0"
WPS, OWS = get_ElementMakerForVersion(VERSION)
xpath_ns = get_xpath_ns(VERSION)

TESTS_HOME = Path(__file__).parent
TD = TESTS_HOME / 'testdata'
CFG_FILE = TESTS_HOME / 'test.cfg'

TESTDATA = {}
TESTDATA['gr4j-cemaneige'] = \
    {'pr': TD / 'gr4j_cemaneige' / 'pr.nc',
     'tas': TD / 'gr4j_cemaneige' / 'tas.nc',
     'evap': TD / 'gr4j_cemaneige' / 'evap.nc'}

TESTDATA['raven-gr4j-cemaneige-nc-ts'] = TD / 'raven-gr4j-cemaneige' / 'Salmon-River-Near-Prince-George_meteo_daily.nc'

TESTDATA['raven-gr4j-cemaneige-nc-rv'] = tuple((TD / 'raven-gr4j-cemaneige').glob('raven-gr4j-salmon.rv?'))

TESTDATA['raven-mohyse-nc-ts'] = TESTDATA['raven-gr4j-cemaneige-nc-ts']
TESTDATA['raven-mohyse'] = TD / 'raven-mohyse'
TESTDATA['raven-mohyse-rv'] = tuple((TD / 'raven-mohyse').glob('raven-mohyse-salmon.rv?'))
TESTDATA['raven-mohyse-ts'] = tuple((TD / 'raven-mohyse').glob('Salmon-River-Near-Prince-George_*.rvt'))

TESTDATA['raven-hmets-nc-ts'] = TESTDATA['raven-gr4j-cemaneige-nc-ts']
TESTDATA['raven-hmets'] = TD / 'raven-hmets'
TESTDATA['raven-hmets-rv'] = tuple((TD / 'raven-hmets').glob('raven-hmets-salmon.rv?'))
TESTDATA['raven-hmets-ts'] = tuple((TD / 'raven-hmets').glob('Salmon-River-Near-Prince-George_*.rvt'))

TESTDATA['raven-hbv-ec-nc-ts'] = TESTDATA['raven-gr4j-cemaneige-nc-ts']
TESTDATA['raven-hbv-ec'] = TD / 'raven-hbv-ec'
TESTDATA['raven-hbv-ec-rv'] = tuple((TD / 'raven-hbv-ec').glob('raven-hbv-ec-salmon.rv?'))
TESTDATA['raven-hbv-ec-ts'] = tuple((TD / 'raven-hbv-ec').glob('Salmon-River-Near-Prince-George_*.rvt'))


class WpsTestClient(WpsClient):

    def get(self, *args, **kwargs):
        query = "?"
        for key, value in kwargs.items():
            query += "{0}={1}&".format(key, value)
        return super(WpsTestClient, self).get(query)


def client_for(service):
    return WpsTestClient(service, WpsTestResponse)


def get_output(doc):
    """Read XML process response and return output dictionary."""
    output = {}
    for output_el in xpath_ns(doc, '/wps:ExecuteResponse'
                                   '/wps:ProcessOutputs/wps:Output'):
        [identifier_el] = xpath_ns(output_el, './ows:Identifier')

        lit_el = xpath_ns(output_el, './wps:Data/wps:LiteralData')
        if lit_el != []:
            output[identifier_el.text] = lit_el[0].text

        ref_el = xpath_ns(output_el, './wps:Reference')
        if ref_el != []:
            output[identifier_el.text] = ref_el[0].attrib['href']

        data_el = xpath_ns(output_el, './wps:Data/wps:ComplexData')
        if data_el != []:
            output[identifier_el.text] = data_el[0].text

    return output


def synthetic_gr4j_inputs(path):
    time = pd.date_range(start='2000-07-01', end='2002-07-01', freq='D')

    pr = 3 * np.ones(len(time))
    pr = xr.DataArray(pr, coords={'time': time}, dims='time', name='pr')
    pr.to_netcdf(os.path.join(path, 'pr.nc'))

    tas = 280 + 20 * np.cos(np.arange(len(time)) * 2 * np.pi / 365.)
    tas = xr.DataArray(tas, coords={'time': time}, dims='time', name='tas')
    tas.to_netcdf(os.path.join(path, 'tas.nc'))

    evap = 3 + 3 * np.cos(-30 + np.arange(len(time)) * 2 * np.pi / 365.)
    evap = xr.DataArray(evap, coords={'time': time}, dims='time', name='evap')
    evap.to_netcdf(os.path.join(path, 'evap.nc'))
