import json

from pywps import Service
from pywps.tests import assert_response_success
from ravenpy.utilities.testdata import get_test_data
from shapely.geometry import MultiPolygon

from raven.processes import NALCMSZonalStatisticsProcess

from .common import CFG_FILE, client_for, count_pixels, get_output


class TestNALCMSZonalStatsProcess:
    def test_simplified_categories(self):
        client = client_for(
            Service(
                processes=[
                    NALCMSZonalStatisticsProcess(),
                ],
                cfgfiles=CFG_FILE,
            )
        )

        fields = [
            "select_all_touching={touches}",
            "simple_categories={simple_categories}",
            "band={band}",
            "shape=file@xlink:href=file://{shape}",
            "raster=file@xlink:href=file://{raster}",
        ]

        datainputs = ";".join(fields).format(
            touches=True,
            simple_categories=True,
            band=1,
            shape=get_test_data("donneesqc_mrc_poly", "mrc_subset.gml")[0],
            raster=get_test_data("cec_nalcms2010_30m", "cec_nalcms_subQC.tiff")[0],
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="nalcms-zonal-stats",
            datainputs=datainputs,
        )

        assert_response_success(resp)
        out = get_output(resp.xml)
        stats = json.loads(out["statistics"])[0]
        assert not {"count", "nodata", "nan"}.issubset(stats)

        geometry = json.loads(out["features"])
        assert isinstance(type(geometry), type(MultiPolygon))

        category_counts = count_pixels(stats)
        assert category_counts == geometry["features"][0]["properties"]["count"]
        assert sum(stats.values()) == geometry["features"][0]["properties"]["count"]

    def test_true_categories(self):
        client = client_for(
            Service(
                processes=[
                    NALCMSZonalStatisticsProcess(),
                ],
                cfgfiles=CFG_FILE,
            )
        )

        fields = [
            "select_all_touching={touches}",
            "simple_categories={simple_categories}",
            "band={band}",
            "shape=file@xlink:href=file://{shape}",
            "raster=file@xlink:href=file://{raster}",
        ]

        datainputs = ";".join(fields).format(
            touches=True,
            simple_categories=False,
            band=1,
            shape=get_test_data("donneesqc_mrc_poly", "mrc_subset.gml")[0],
            raster=get_test_data("cec_nalcms2010_30m", "cec_nalcms_subQC.tiff")[0],
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="nalcms-zonal-stats",
            datainputs=datainputs,
        )

        assert_response_success(resp)
        out = get_output(resp.xml)
        stats = json.loads(out["statistics"])[0]
        assert not {"count", "nodata", "nan"}.issubset(stats)

        geometry = json.loads(out["features"])
        assert isinstance(type(geometry), type(MultiPolygon))

        category_counts = count_pixels(stats)
        assert category_counts == geometry["features"][0]["properties"]["count"]
        assert sum(stats.values()) == geometry["features"][0]["properties"]["count"]

    def test_wcs_simplified_categories(self):
        client = client_for(
            Service(
                processes=[
                    NALCMSZonalStatisticsProcess(),
                ],
                cfgfiles=CFG_FILE,
            )
        )
        fields = [
            "select_all_touching={touches}",
            "simple_categories={simple_categories}",
            "band={band}",
            "shape=file@xlink:href=file://{shape}",
        ]

        datainputs = ";".join(fields).format(
            touches=True,
            simple_categories=True,
            band=1,
            shape=get_test_data("watershed_vector", "Basin_test.zip")[0],
        )
        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="nalcms-zonal-stats",
            datainputs=datainputs,
        )

        assert_response_success(resp)
        out = get_output(resp.xml)
        stats = json.loads(out["statistics"])[0]
        assert not {"count", "nodata", "nan"}.issubset(stats)

        geometry = json.loads(out["features"])
        assert isinstance(type(geometry), type(MultiPolygon))

        category_counts = count_pixels(stats)
        assert category_counts == geometry["features"][0]["properties"]["count"]
        assert sum(stats.values()) == geometry["features"][0]["properties"]["count"]

    def test_wcs_true_categories(self):
        client = client_for(
            Service(
                processes=[
                    NALCMSZonalStatisticsProcess(),
                ],
                cfgfiles=CFG_FILE,
            )
        )
        fields = [
            "select_all_touching={touches}",
            "simple_categories={simple_categories}",
            "band={band}",
            "shape=file@xlink:href=file://{shape}",
        ]

        datainputs = ";".join(fields).format(
            touches=True,
            simple_categories=False,
            band=1,
            shape=get_test_data("watershed_vector", "Basin_test.zip")[0],
        )
        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="nalcms-zonal-stats",
            datainputs=datainputs,
        )

        assert_response_success(resp)
        out = get_output(resp.xml)
        stats = json.loads(out["statistics"])[0]
        assert not {"count", "nodata", "nan"}.issubset(stats)

        geometry = json.loads(out["features"])
        assert isinstance(type(geometry), type(MultiPolygon))

        category_counts = count_pixels(stats)
        assert category_counts == geometry["features"][0]["properties"]["count"]
        assert sum(stats.values()) == geometry["features"][0]["properties"]["count"]
