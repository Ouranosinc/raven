from pywps import Service
from pywps.tests import assert_response_success
from .common import client_for, TESTDATA, CFG_FILE, get_output, urlretrieve

from raven.processes import ShapeAreaProcess


class TestShapeAreaProcess:

    def test_simple(self):
        client = client_for(Service(processes=[ShapeAreaProcess(), ], cfgfiles=CFG_FILE))

        fields = ['use_all_features={feat}', 'crs={crs}', 'shape=file@xlink:href=file://{file}']
        datainputs = ';'.join(fields).format(
            feat=True,
            crs='4326',
            file=TESTDATA['watershed_vector']
        )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='shape-area', datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)

        assert 'area' in out
        assert 'centroids' in out
        assert 'schemas' in out
