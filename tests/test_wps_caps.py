from pywps import Service
from raven.processes import processes

from .common import client_for


def test_wps_caps():
    client = client_for(Service(processes=processes))
    resp = client.get(service="wps", request="getcapabilities", version="1.0.0")
    names = resp.xpath_text(
        "/wps:Capabilities" "/wps:ProcessOfferings" "/wps:Process" "/ows:Identifier"
    )
    sn = set(names.split())

    assert sn == {
        "gr4j-cemaneige",
        "raster-subset",
        "raven",
        "raven-multi-model",
        "raven-gr4j-cemaneige",
        "raven-mohyse",
        "raven-hmets",
        "raven-hbv-ec",
        "shape-properties",
        "hydrobasins-select",
        "terrain-analysis",
        "nalcms-zonal-stats",
        "nalcms-zonal-stats-raster",
        "zonal-stats",
        "ostrich-gr4j-cemaneige",
        "ostrich-mohyse",
        "ostrich-hmets",
        "ostrich-hbv-ec",
        "objective-function",
        "graph_ensemble_uncertainty",
        "graph_single_hydrograph",
        "ts_stats",
        "freq_analysis",
        "base_flow_index",
        "ts_stats_graph",
        "regionalisation",
        "hindcast-evaluation",
        "graph_objective_function_fit",
        "fit",
        "graph_fit",
        "climatology_esp",
        "forecast-floodrisk",
        "hindcasting",
        "realtime-forecast",
    }
