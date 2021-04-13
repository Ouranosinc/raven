# -*- coding: utf-8 -*-
import json
import logging
import tempfile

import geopandas as gpd
from pywps import FORMATS, ComplexOutput, LiteralInput, Process
from pywps.exceptions import InvalidParameterValue
from ravenpy.utilities import geoserver
from ravenpy.utilities.checks import single_file_check
from ravenpy.utilities.io import archive_sniffer, crs_sniffer

from raven.utils import parse_lonlat

LOGGER = logging.getLogger("PYWPS")


class HydroBasinsSelectionProcess(Process):
    """Given lat/lon coordinates that point to a North American watershed,
    return the feature containing the coordinates or the entire upstream water basin."""

    def __init__(self):
        inputs = [
            LiteralInput(
                "location",
                "Location coordinates (lon, lat)",
                abstract="Location coordinates (longitude, latitude) for point of interest.",
                data_type="string",
                min_occurs=1,
                max_occurs=2,
            ),
            # LiteralInput('level', 'Resolution level of HydroBASINS Shapes',
            #              data_type='integer',
            #              default=12,
            #              allowed_values=[7, 8, 9, 10, 11, 12],
            #              min_occurs=0,
            #              max_occurs=1),
            # LiteralInput('lakes', 'Use the HydroBASINS version that includes lake outlines',
            #              data_type='boolean',
            #              default='true',
            #              min_occurs=0,
            #              max_occurs=1),
            LiteralInput(
                "aggregate_upstream",
                "Attempt to capture both the containing basin and all tributary basins from point",
                data_type="boolean",
                default="false",
                min_occurs=0,
                max_occurs=1,
            ),
        ]

        outputs = [
            ComplexOutput(
                "feature",
                "Watershed feature geometry",
                abstract="Geographic representation of shape properties.",
                supported_formats=[FORMATS.GEOJSON],
            ),
            ComplexOutput(
                "upstream_ids",
                "HydroBASINS IDs for all immediate upstream basins",
                abstract="List of all tributary sub-basins according to their HydroBASINS IDs, "
                "including the downstream basin.",
                supported_formats=[FORMATS.JSON],
            ),
        ]

        super(HydroBasinsSelectionProcess, self).__init__(
            self._handler,
            identifier="hydrobasins-select",
            title="Select a HydroBASINS watershed geometry",
            version="1.1",
            abstract="Return a watershed from the HydroSheds database as a polygon vector file.",
            metadata=[],
            inputs=inputs,
            outputs=outputs,
            status_supported=True,
            store_supported=True,
        )

    def _handler(self, request, response):

        # level = 12  # request.inputs['level'][0].data
        # lakes = True  # request.inputs['lakes'][0].data
        collect_upstream = request.inputs["aggregate_upstream"][0].data
        lon, lat = parse_lonlat(request.inputs["location"][0].data)

        bbox = (lon, lat, lon, lat)

        domain = geoserver.select_hybas_domain(bbox)
        hybas_request = geoserver.get_hydrobasins_location_wfs(bbox, domain=domain)

        response.update_status("Found downstream watershed", status_percentage=20)

        # Find HYBAS_ID
        gdf = gpd.read_file(hybas_request.decode())
        id_number = gdf["id"][0]

        if collect_upstream:

            # Collect features from GeoServer
            response.update_status("Collecting relevant features", status_percentage=50)

            # Identify upstream sub-basins and write to a new file
            upstream = geoserver.hydrobasins_upstream(gdf.loc[0], domain)

            # Aggregate upstream features into a single geometry.
            agg = geoserver.hydrobasins_aggregate(upstream)

            response.update_status("Writing upstream features", status_percentage=70)
            # The aggregation returns a FeatureCollection with one feature. We select the first feature so that the
            # output is a Feature whether aggregate is True or False.
            afeat = json.loads(agg.to_json())["features"][0]
            response.outputs["feature"].data = json.dumps(afeat)
            response.outputs["upstream_ids"].data = json.dumps(upstream["id"].tolist())

        else:
            response.outputs["feature"].data = gdf.loc[0].to_json()
            response.outputs["upstream_ids"].data = json.dumps([id_number])

        return response
