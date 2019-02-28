from pywps import Service
from pywps.tests import assert_response_success, assert_process_exception
from .common import client_for, TESTDATA, CFG_FILE, get_output
import json
from raven.processes import ShapePropertiesProcess


class TestGenericShapePropertiesProcess:

    def test_simple(self):
        client = client_for(Service(processes=[ShapePropertiesProcess(), ], cfgfiles=CFG_FILE))

        fields = ['shape=file@xlink:href=file://{file}', 'crs={crs}', 'projected_crs={projected_crs}']
        datainputs = ';'.join(fields).format(
            file=TESTDATA['watershed_vector'],
            crs=4326,
            projected_crs=32198
        )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='shape-properties', datainputs=datainputs)

        assert_response_success(resp)

        out = get_output(resp.xml)
        assert {'properties'}.issubset([*out])

        props = json.loads(out['properties'])
        assert {'centroid', 'area', 'perimeter', 'gravelius'}.issubset(props[0].keys())

    def test_bad_epsg(self):
        client = client_for(Service(processes=[ShapePropertiesProcess(), ], cfgfiles=CFG_FILE))

        fields = ['shape=file@xlink:href=file://{file}', 'crs={crs}', 'projected_crs={projected_crs}']
        datainputs = ';'.join(fields).format(
            file=TESTDATA['watershed_vector'],
            crs=4326,
            projected_crs=4326
        )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='shape-properties', datainputs=datainputs)

        assert_response_success(resp)

        out = get_output(resp.xml)
        assert {'properties'}.issubset([*out])

        props = json.loads(out['properties'])
        assert {'centroid', 'area', 'perimeter', 'gravelius'}.issubset(props[0].keys())
