# -*- coding: utf-8 -*-
import json
import logging
import tempfile

import fiona
import geopandas as gpd
from pywps import LiteralInput, ComplexOutput
from pywps import Process, FORMATS
from pywps.exceptions import InvalidParameterValue

from raven.utilities import gis
from raven.utils import archive_sniffer, single_file_check, parse_lonlat
from raven.utils import crs_sniffer

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
            LiteralInput(
                "method",
                "WFS request method to use for passing data.",
                data_type="string",
                default="GET",
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
            version="1.2",
            abstract="Return a watershed from the HydroSheds database as a polygon vector file.",
            metadata=[],
            inputs=inputs,
            outputs=outputs,
            status_supported=True,
            store_supported=True,
        )

    def _handler(self, request, response):

        level = 12  # request.inputs['level'][0].data
        lakes = True  # request.inputs['lakes'][0].data
        collect_upstream = request.inputs["aggregate_upstream"][0].data
        lon, lat = parse_lonlat(request.inputs["location"][0].data)
        method = request.inputs["method"][0].data

        bbox = (lon, lat, lon, lat)

        shape_url = tempfile.NamedTemporaryFile(
            prefix="hybas_", suffix=".gml", delete=False, dir=self.workdir
        ).name

        domain = gis.select_hybas_domain(bbox)
        hybas_gml = gis.get_hydrobasins_location_wfs(
            coordinates=bbox, lakes=lakes, level=level, domain=domain, method=method
        )

        if isinstance(hybas_gml, str):
            write_flags = "w"
        else:
            write_flags = "wb"

        with open(shape_url, write_flags) as f:
            f.write(hybas_gml)

        response.update_status("Found downstream watershed", status_percentage=10)

        extensions = [".gml", ".shp", ".gpkg", ".geojson", ".json"]
        shp = single_file_check(
            archive_sniffer(shape_url, working_dir=self.workdir, extensions=extensions)
        )

        shape_crs = crs_sniffer(shp)

        # Find HYBAS_ID
        src = fiona.open(shp, "r", crs=shape_crs)
        feat = next(iter(src))
        hybas_id = feat["properties"]["HYBAS_ID"]
        gml_id = feat["properties"]["gml_id"]

        if collect_upstream:

            main_bas = feat["properties"]["MAIN_BAS"]

            if lakes is False or level != 12:
                raise InvalidParameterValue("Set lakes to True and level to 12.")

            # Collect features from GeoServer
            response.update_status("Collecting relevant features", status_percentage=70)

            region_url = gis.get_hydrobasins_attributes_wfs(
                attribute="MAIN_BAS",
                value=main_bas,
                lakes=lakes,
                level=level,
                domain=domain,
            )

            # Read table of relevant features sharing main basin
            df = gpd.read_file(region_url)

            # TODO: Load and keep this data in memory; Figure out how to better handle encoding and column names.
            # Identify upstream sub-basins and write to a new file
            up = gis.hydrobasins_upstream_ids(hybas_id, df)
            upfile = tempfile.NamedTemporaryFile(
                prefix="hybas_", suffix=".json", delete=False, dir=self.workdir
            ).name
            up.to_file(upfile, driver="GeoJSON")

            # Aggregate upstream features into a single geometry.
            gdf = gpd.read_file(upfile)
            agg = gis.hydrobasins_aggregate(gdf)

            # The aggregation returns a FeatureCollection with one feature. We select the first feature so that the
            # output is a Feature whether aggregate is True or False.
            afeat = json.loads(agg.to_json())["features"][0]
            response.outputs["feature"].data = json.dumps(afeat)
            response.outputs["upstream_ids"].data = json.dumps(up["id"].tolist())

        else:
            response.outputs["feature"].data = json.dumps(feat)
            response.outputs["upstream_ids"].data = json.dumps([gml_id,])

        src.close()

        return response
