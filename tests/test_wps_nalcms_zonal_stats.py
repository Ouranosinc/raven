import json

import pytest
from common import CFG_FILE, client_for, count_pixels, get_output
from metalink import download as md
from pywps import Service, get_ElementMakerForVersion
from pywps.app.basic import get_xpath_ns
from pywps.tests import assert_response_success
from shapely.geometry import MultiPolygon

from raven.processes import (
    NALCMSZonalStatisticsProcess,
    NALCMSZonalStatisticsRasterProcess,
)

VERSION = "1.0.0"
WPS, OWS = get_ElementMakerForVersion(VERSION)
xpath_ns = get_xpath_ns(VERSION)


class TestNALCMSZonalStatsProcess:
    def test_simplified_categories(self, yangtze):
        client = client_for(
            Service(processes=[NALCMSZonalStatisticsProcess()], cfgfiles=CFG_FILE)
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
            shape=yangtze.fetch("donneesqc_mrc_poly/mrc_subset.gml"),
            raster=yangtze.fetch("cec_nalcms2010_30m/cec_nalcms_subQC.tiff"),
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

    def test_true_categories(self, yangtze):
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
            shape=yangtze.fetch("donneesqc_mrc_poly/mrc_subset.gml"),
            raster=yangtze.fetch("cec_nalcms2010_30m/cec_nalcms_subQC.tiff"),
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

    def test_wcs_simplified_categories(self, yangtze):
        client = client_for(
            Service(processes=[NALCMSZonalStatisticsProcess()], cfgfiles=CFG_FILE)
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
            shape=yangtze.fetch("watershed_vector/Basin_test.zip"),
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

    def test_wcs_true_categories(self, yangtze):
        client = client_for(
            Service(processes=[NALCMSZonalStatisticsProcess()], cfgfiles=CFG_FILE)
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
            shape=yangtze.fetch("watershed_vector/Basin_test.zip"),
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


@pytest.mark.online
class TestNALCMSZonalStatsWithRasterProcess:
    def test_wcs_simplified_categories(self, yangtze):
        client = client_for(
            Service(processes=[NALCMSZonalStatisticsRasterProcess()], cfgfiles=CFG_FILE)
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
            shape=yangtze.fetch("watershed_vector/Basin_test.zip"),
        )
        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="nalcms-zonal-stats-raster",
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

        assert {"raster"}.issubset([*out])
        d = md.get(out["raster"], path="/tmp", segmented=False)
        assert {f"/tmp/subset_{n}.tiff" for n in range(1, 6)}.issubset(set(d))

    def test_mississippi(self, yangtze):
        """Extremely large region testing timeout value for the GeoServer response."""
        client = client_for(
            Service(processes=[NALCMSZonalStatisticsRasterProcess()], cfgfiles=CFG_FILE)
        )
        fields = [
            "shape=file@xlink:href=file://{shape}",
        ]

        datainputs = ";".join(fields).format(
            shape=yangtze.fetch("polygons/mississippi.geojson")
        )
        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="nalcms-zonal-stats-raster",
            datainputs=datainputs,
        )
        el = resp.xpath("//wps:ExceptionReport/ows:Exception/ows:ExceptionText")[0]
        assert "areas smaller" in el.text
