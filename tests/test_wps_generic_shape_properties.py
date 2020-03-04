import geojson
import json
import tempfile
import numpy as np
from pywps import Service
from pywps.tests import assert_response_success

from raven.processes import ShapePropertiesProcess
from .common import client_for, TESTDATA, CFG_FILE, get_output


class TestGenericShapePropertiesProcess:

    @staticmethod
    def make_shape():
        coords = [
            [
                [
                    -74.827880859375,
                    46.240651955001695
                ],
                [
                    -74.3280029296875,
                    46.027481852486645
                ],
                [
                    -73.7786865234375,
                    46.02366774426006
                ],
                [
                    -73.1195068359375,
                    46.12274903582433
                ],
                [
                    -72.7130126953125,
                    46.32796494040746
                ],
                [
                    -72.6690673828125,
                    46.694667307773116
                ],
                [
                    -73.037109375,
                    46.86770273172814
                ],
                [
                    -73.5589599609375,
                    46.87145819560722
                ],
                [
                    -73.927001953125,
                    46.78877728793222
                ],
                [
                    -74.080810546875,
                    46.56641407568593
                ],
                [
                    -73.86657714843749,
                    46.437856895024204
                ],
                [
                    -73.5479736328125,
                    46.47191632087041
                ],
                [
                    -73.267822265625,
                    46.58906908309182
                ],
                [
                    -73.6798095703125,
                    46.649436163350245
                ],
                [
                    -73.333740234375,
                    46.717268685073954
                ],
                [
                    -72.9876708984375,
                    46.55130547880643
                ],
                [
                    -73.267822265625,
                    46.392411189814645
                ],
                [
                    -73.62487792968749,
                    46.1912395780416
                ],
                [
                    -74.1357421875,
                    46.0998999106273
                ],
                [
                    -74.827880859375,
                    46.240651955001695
                ]
            ]
        ]

        geo_def = geojson.Polygon(coords)
        # raise Warning(geo_def)
        # raise Exception(geo_def.errors())

        temp = tempfile.NamedTemporaryFile(suffix='.json', delete=False)
        with open(temp.name, 'w') as f:
            geojson.dump(geo_def, f, indent=4)
        return temp

    def test_simple(self):
        shape = self.make_shape()
        client = client_for(Service(processes=[ShapePropertiesProcess(), ], cfgfiles=CFG_FILE))

        fields = ['shape=file@xlink:href=file://{file}', 'crs={crs}', 'projected_crs={projected_crs}']
        datainputs = ';'.join(fields).format(
            file=shape.name,
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

        # np.testing.assert_approx_equal(props[0]['area'], 24804229364.)
        # np.testing.assert_approx_equal(props[0]['centroid'][0], -73.58396746)
        # np.testing.assert_approx_equal(props[0]['centroid'][1], 46.88002938)
        # np.testing.assert_approx_equal(props[0]['perimeter'], 1520538.80767662)
        # np.testing.assert_approx_equal(props[0]['gravelius'], 2.7235145984)

    def test_geographic_epsg(self):
        """Calculate the geometric properties using degree-length units"""

        shape = self.make_shape()
        client = client_for(Service(processes=[ShapePropertiesProcess(), ], cfgfiles=CFG_FILE))

        fields = ['shape=file@xlink:href=file://{file}', 'crs={crs}', 'projected_crs={projected_crs}']
        datainputs = ';'.join(fields).format(
            file=shape.name,
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

        # np.testing.assert_approx_equal(props[0]['area'], 2.9369078)
        # np.testing.assert_approx_equal(props[0]['centroid'][0], -73.58396746)
        # np.testing.assert_approx_equal(props[0]['centroid'][1], 46.88002938)
        # np.testing.assert_approx_equal(props[0]['perimeter'], 18.405642329)
        # np.testing.assert_approx_equal(props[0]['gravelius'], 3.02970880)

    def test_multifeature_geojson(self):
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

    def test_multifeature_zipped_shapefile(self):
        """Calculate shape properties for multiple features in a shape"""

        client = client_for(Service(processes=[ShapePropertiesProcess(), ], cfgfiles=CFG_FILE))

        fields = ['shape=file@xlink:href=file://{file}', 'crs={crs}', 'projected_crs={projected_crs}']
        datainputs = ';'.join(fields).format(
            file=TESTDATA['mrc_subset_zipped'],
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
