import json
import numpy as np

from pywps import Service
from pywps.tests import assert_response_success
from .common import client_for, TESTDATA, CFG_FILE, get_output
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

        np.testing.assert_approx_equal(props[0]['area'], 44877188052)
        np.testing.assert_approx_equal(props[0]['centroid'][0], -72.69128332)
        np.testing.assert_approx_equal(props[0]['centroid'][1], 49.50119363)
        np.testing.assert_approx_equal(props[0]['perimeter'], 1861151.69970881)
        np.testing.assert_approx_equal(props[0]['gravelius'], 2.4783578)

    def test_geographic_epsg(self):
        """Calculate the geometric properties using degree-length units"""

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

        np.testing.assert_approx_equal(props[0]['area'], 5.63123920)
        np.testing.assert_approx_equal(props[0]['centroid'][0], -72.69128332)
        np.testing.assert_approx_equal(props[0]['centroid'][1], 49.50119363)
        np.testing.assert_approx_equal(props[0]['perimeter'], 20.464501)
        np.testing.assert_approx_equal(props[0]['gravelius'], 2.4327318)

    def test_multifeature_shape(self):
        """Calculate shape properties for multiple features in a shape"""

        client = client_for(Service(processes=[ShapePropertiesProcess(), ], cfgfiles=CFG_FILE))

        fields = ['shape=file@xlink:href=file://{file}', 'crs={crs}', 'projected_crs={projected_crs}']
        datainputs = ';'.join(fields).format(
            file=TESTDATA['mrc_subset'],
            crs=4326,
            projected_crs=32198
        )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='shape-properties', datainputs=datainputs)

        assert_response_success(resp)

        out = get_output(resp.xml)
        assert {'properties'}.issubset([*out])

        props = json.loads(out['properties'])
        for i in range(len(props)):
            assert {'centroid', 'area', 'perimeter', 'gravelius'}.issubset(props[i].keys())

        np.testing.assert_approx_equal(props[0]['area'], 111417901.6141605)
        np.testing.assert_approx_equal(props[0]['centroid'][0], -71.8223648)
        np.testing.assert_approx_equal(props[0]['centroid'][1], 48.8974365)
        np.testing.assert_approx_equal(props[0]['perimeter'], 46351.1628725)
        np.testing.assert_approx_equal(props[0]['gravelius'], 1.2387344)

        np.testing.assert_approx_equal(props[-1]['area'], 334136220.2845515)
        np.testing.assert_approx_equal(props[-1]['centroid'][0], -72.6117018)
        np.testing.assert_approx_equal(props[-1]['centroid'][1], 46.3632907)
        np.testing.assert_approx_equal(props[-1]['perimeter'], 92477.3057504)
        np.testing.assert_approx_equal(props[-1]['gravelius'], 1.4271461)
