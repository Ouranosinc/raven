import json
import tempfile

import geojson
import numpy as np
from common import CFG_FILE, client_for, get_output
from pywps import Service
from pywps.tests import assert_response_success

from raven.processes import ShapePropertiesProcess


class TestGenericShapePropertiesProcess:
    @staticmethod
    def make_shape():
        coords = [
            [
                [-74.827880859375, 46.240651955001695],
                [-74.3280029296875, 46.027481852486645],
                [-73.7786865234375, 46.02366774426006],
                [-73.1195068359375, 46.12274903582433],
                [-72.7130126953125, 46.32796494040746],
                [-72.6690673828125, 46.694667307773116],
                [-73.037109375, 46.86770273172814],
                [-73.5589599609375, 46.87145819560722],
                [-73.927001953125, 46.78877728793222],
                [-74.080810546875, 46.56641407568593],
                [-73.86657714843749, 46.437856895024204],
                [-73.5479736328125, 46.47191632087041],
                [-73.267822265625, 46.58906908309182],
                [-73.6798095703125, 46.649436163350245],
                [-73.333740234375, 46.717268685073954],
                [-72.9876708984375, 46.55130547880643],
                [-73.267822265625, 46.392411189814645],
                [-73.62487792968749, 46.1912395780416],
                [-74.1357421875, 46.0998999106273],
                [-74.827880859375, 46.240651955001695],
            ]
        ]

        geo_def = geojson.Polygon(coords)
        # raise Warning(geo_def)
        # raise Exception(geo_def.errors())

        temp = tempfile.NamedTemporaryFile(suffix=".json", delete=False)
        with open(temp.name, "w") as f:
            geojson.dump(geo_def, f, indent=4)
        return temp

    def test_simple(self):
        shape = self.make_shape()
        client = client_for(
            Service(processes=[ShapePropertiesProcess()], cfgfiles=CFG_FILE)
        )

        fields = [
            "shape=file@xlink:href=file://{file}",
            "crs={crs}",
            "projected_crs={projected_crs}",
        ]
        datainputs = ";".join(fields).format(
            file=shape.name, crs=4326, projected_crs=6622
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="shape-properties",
            datainputs=datainputs,
        )

        assert_response_success(resp)

        out = get_output(resp.xml)
        assert {"properties"}.issubset([*out])

        props = json.loads(out["properties"])
        assert {"centroid", "area", "perimeter", "gravelius"}.issubset(props[0].keys())

        np.testing.assert_allclose(props[0]["perimeter"], 673431, atol=1)
        np.testing.assert_approx_equal(props[0]["area"], 6258366698.5253, 4)
        np.testing.assert_approx_equal(props[0]["centroid"][0], -73.41117680)
        np.testing.assert_approx_equal(props[0]["centroid"][1], 46.46286765)
        np.testing.assert_approx_equal(props[0]["gravelius"], 2.4013618703)

    def test_geographic_epsg(self):
        """Calculate the geometric properties using degree-length units"""

        shape = self.make_shape()
        client = client_for(
            Service(processes=[ShapePropertiesProcess()], cfgfiles=CFG_FILE)
        )

        fields = [
            "shape=file@xlink:href=file://{file}",
            "crs={crs}",
            "projected_crs={projected_crs}",
        ]
        datainputs = ";".join(fields).format(
            file=shape.name, crs=4326, projected_crs=4326
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="shape-properties",
            datainputs=datainputs,
        )

        assert_response_success(resp)

        out = get_output(resp.xml)
        assert {"properties"}.issubset([*out])

        props = json.loads(out["properties"])
        assert {"centroid", "area", "perimeter", "gravelius"}.issubset(props[0].keys())

        np.testing.assert_approx_equal(props[0]["area"], 0.7342578)
        np.testing.assert_approx_equal(props[0]["centroid"][0], -73.41117680)
        np.testing.assert_approx_equal(props[0]["centroid"][1], 46.46286765)
        np.testing.assert_approx_equal(props[0]["perimeter"], 8.14412255)
        np.testing.assert_approx_equal(props[0]["gravelius"], 2.68111268)

    def test_multifeature_geojson(self, yangtze):
        """Calculate shape properties for multiple features in a shape"""

        client = client_for(
            Service(processes=[ShapePropertiesProcess()], cfgfiles=CFG_FILE)
        )

        fields = [
            "shape=file@xlink:href=file://{file}",
            "crs={crs}",
            "projected_crs={projected_crs}",
        ]
        datainputs = ";".join(fields).format(
            file=yangtze.fetch("donneesqc_mrc_poly/mrc_subset.gml"),
            crs=4326,
            projected_crs=6622,
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="shape-properties",
            datainputs=datainputs,
        )

        assert_response_success(resp)

        out = get_output(resp.xml)
        assert {"properties"}.issubset([*out])

        props = json.loads(out["properties"])
        for i in range(len(props)):
            assert {"centroid", "area", "perimeter", "gravelius"}.issubset(
                props[i].keys()
            )

        np.testing.assert_allclose(props[0]["area"], 111417901, atol=1)
        np.testing.assert_approx_equal(props[0]["centroid"][0], -71.8223648)
        np.testing.assert_approx_equal(props[0]["centroid"][1], 48.8974365)
        np.testing.assert_approx_equal(props[0]["gravelius"], 1.2387344)
        np.testing.assert_approx_equal(props[0]["perimeter"], 46351.1628725)

        np.testing.assert_allclose(props[-1]["area"], 334136220, atol=100)
        np.testing.assert_allclose(props[-1]["perimeter"], 92477.3, atol=0.1)
        np.testing.assert_approx_equal(props[-1]["centroid"][0], -72.6117018)
        np.testing.assert_approx_equal(props[-1]["centroid"][1], 46.3632907)
        np.testing.assert_approx_equal(props[-1]["gravelius"], 1.4271461)

    def test_multifeature_zipped_shapefile(self, yangtze):
        """Calculate shape properties for multiple features in a shape"""

        client = client_for(
            Service(processes=[ShapePropertiesProcess()], cfgfiles=CFG_FILE)
        )

        fields = [
            "shape=file@xlink:href=file://{file}",
            "crs={crs}",
            "projected_crs={projected_crs}",
        ]
        datainputs = ";".join(fields).format(
            file=yangtze.fetch("donneesqc_mrc_poly/mrc_subset.zip"),
            crs=4326,
            projected_crs=6622,
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="shape-properties",
            datainputs=datainputs,
        )

        assert_response_success(resp)

        out = get_output(resp.xml)
        assert {"properties"}.issubset([*out])

        props = json.loads(out["properties"])
        for i in range(len(props)):
            assert {"centroid", "area", "perimeter", "gravelius"}.issubset(
                props[i].keys()
            )

        np.testing.assert_allclose(props[0]["area"], 111417901, atol=1)
        np.testing.assert_approx_equal(props[0]["centroid"][0], -71.8223648)
        np.testing.assert_approx_equal(props[0]["centroid"][1], 48.8974365)
        np.testing.assert_approx_equal(props[0]["gravelius"], 1.2387344)
        np.testing.assert_approx_equal(props[0]["perimeter"], 46351.1628725)

        np.testing.assert_allclose(props[-1]["area"], 334136220, atol=1)
        np.testing.assert_allclose(props[-1]["perimeter"], 92477.3, atol=0.1)
        np.testing.assert_approx_equal(props[-1]["centroid"][0], -72.6117018)
        np.testing.assert_approx_equal(props[-1]["centroid"][1], 46.3632907)
        np.testing.assert_approx_equal(props[-1]["gravelius"], 1.4271461)
