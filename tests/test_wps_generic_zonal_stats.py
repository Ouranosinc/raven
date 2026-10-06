import json

import pytest
from common import CFG_FILE, client_for, count_pixels, get_output
from pywps import Service
from pywps.tests import assert_response_success
from shapely.geometry import MultiPolygon, shape

from raven.processes import ZonalStatisticsProcess


class TestGenericZonalStatsProcess:
    def test_simple(self, yangtze):
        client = client_for(
            Service(
                processes=[
                    ZonalStatisticsProcess(),
                ],
                cfgfiles=CFG_FILE,
            )
        )

        fields = [
            "select_all_touching=True",
            "categorical=False",
            "band=1",
            f"shape=file@xlink:href=file://{yangtze.fetch('donneesqc_mrc_poly/mrc_subset.gml')}",
            f"raster=file@xlink:href=file://{yangtze.fetch('earthenv_dem_90m/earthenv_dem90_southernQuebec.tiff')}",
        ]
        datainputs = ";".join(fields)

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="zonal-stats",
            datainputs=datainputs,
        )

        assert_response_success(resp)
        out = get_output(resp.xml)
        feature = json.loads(out["statistics"])["features"][0]

        stats = feature["properties"]
        assert {"count", "min", "max", "mean", "median", "sum", "nodata"}.issubset(
            stats
        )

        geometry = shape(feature["geometry"])
        assert isinstance(type(geometry), type(MultiPolygon))

    def test_categorized(self, yangtze):
        client = client_for(
            Service(
                processes=[
                    ZonalStatisticsProcess(),
                ],
                cfgfiles=CFG_FILE,
            )
        )

        fields = [
            "select_all_touching=True",
            "categorical=True",
            "band=1",
            f"shape=file@xlink:href=file://{yangtze.fetch('donneesqc_mrc_poly/mrc_subset.gml')}",
            "raster=file@xlink:href=file://{yangtze.fetch('earthenv_dem_90m/earthenv_dem90_southernQuebec.tiff')",
        ]
        datainputs = ";".join(fields)

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="zonal-stats",
            datainputs=datainputs,
        )

        assert_response_success(resp)
        out = get_output(resp.xml)
        feature = json.loads(out["statistics"])["features"][0]

        stats = feature["properties"]
        assert {"count", "min", "max", "mean", "median", "sum", "nodata"}.issubset(
            stats
        )

        # Check for accurate pixel counts
        category_counts = count_pixels(stats, numeric_categories=True)
        assert category_counts == stats["count"]

        geometry = shape(feature["geometry"])
        assert isinstance(type(geometry), type(MultiPolygon))

    @pytest.mark.online
    @pytest.mark.xfail(
        raises=AssertionError, reason="Geoserver unavailable", strict=False
    )
    def test_geoserver_dem_wcs_simple(self, yangtze):
        client = client_for(
            Service(
                processes=[
                    ZonalStatisticsProcess(),
                ],
                cfgfiles=CFG_FILE,
            )
        )
        fields = [
            "select_all_touching=True",
            "categorical=False",
            "band=1",
            f"shape=file@xlink:href=file://{yangtze.fetch('donneesqc_mrc_poly/mrc_subset.gml')}",
        ]
        datainputs = ";".join(fields)

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="zonal-stats",
            datainputs=datainputs,
        )

        assert_response_success(resp)
        out = get_output(resp.xml)
        feature = json.loads(out["statistics"])["features"][0]

        stats = feature["properties"]
        assert {"count", "min", "max", "mean", "median", "sum", "nodata"}.issubset(
            stats
        )

        geometry = shape(feature["geometry"])
        assert isinstance(type(geometry), type(MultiPolygon))

    @pytest.mark.online
    @pytest.mark.xfail(
        raises=AssertionError, reason="Geoserver unavailable", strict=False
    )
    def test_geoserver_dem_wcs_categorized(self, yangtze):
        client = client_for(
            Service(
                processes=[
                    ZonalStatisticsProcess(),
                ],
                cfgfiles=CFG_FILE,
            )
        )
        fields = [
            "select_all_touching=True",
            "categorical=True",
            "band=1",
            f"shape=file@xlink:href=file://{yangtze.fetch('donneesqc_mrc_poly/mrc_subset.gml')}",
        ]
        datainputs = ";".join(fields)

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="zonal-stats",
            datainputs=datainputs,
        )

        assert_response_success(resp)
        out = get_output(resp.xml)
        feature = json.loads(out["statistics"])["features"][0]

        stats = feature["properties"]
        assert {"count", "min", "max", "mean", "median", "sum", "nodata"}.issubset(
            stats
        )

        # Check for accurate pixel counts
        category_counts = count_pixels(stats, numeric_categories=True)
        assert category_counts == stats["count"]

        geometry = shape(feature["geometry"])
        assert isinstance(type(geometry), type(MultiPolygon))
