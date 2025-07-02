import json
import logging

import fiona
from pyproj.crs import CRS
from pywps import FORMATS, ComplexInput, ComplexOutput, LiteralInput, Process
from shapely.geometry import shape

from raven.utilities.analysis import geom_prop
from raven.utilities.checks import multipolygon_check, single_file_check
from raven.utilities.geo import geom_transform
from raven.utilities.io import archive_sniffer, crs_sniffer

LOGGER = logging.getLogger("PYWPS")


class ShapePropertiesProcess(Process):
    """Given a file containing vector data, provide general information and spatial characteristics."""

    def __init__(self):
        inputs = [
            ComplexInput(
                "shape",
                "Vector Shape",
                abstract="An ESRI Shapefile, GML, GeoPackage, JSON or GeoJSON file."
                " The ESRI Shapefile must be zipped and contain the .shp, .shx, and .dbf.",
                supported_formats=[
                    FORMATS.GML,
                    FORMATS.GEOJSON,
                    FORMATS.SHP,
                    FORMATS.JSON,
                    FORMATS.ZIP,
                ],
                min_occurs=1,
                max_occurs=1,
            ),
            LiteralInput(
                "projected_crs",
                "Coordinate Reference System for area calculation (Default: EPSG:6622, NAD83(CSRS) / Quebec Lambert)",
                data_type="integer",
                default=6622,
                min_occurs=1,
                max_occurs=1,
            ),
        ]
        outputs = [
            ComplexOutput(
                "properties",
                "Feature schemas",
                abstract="Geographic representations and descriptions of shape properties: "
                "centroid coordinates, area, perimeter and gravelius index.",
                supported_formats=[FORMATS.JSON],
            ),
        ]

        super().__init__(
            self._handler,
            identifier="shape-properties",
            title="Shape Properties",
            version="1.0",
            abstract="Return shape area in square metres based on line boundaries of a polygonal vector file.",
            metadata=[],
            inputs=inputs,
            outputs=outputs,
            status_supported=True,
            store_supported=True,
        )

    def _handler(self, request, response):
        shape_url = request.inputs["shape"][0].file
        projected_crs = request.inputs["projected_crs"][0].data

        extensions = [".gml", ".shp", ".gpkg", ".geojson", ".json"]
        vector_file = single_file_check(
            archive_sniffer(shape_url, working_dir=self.workdir, extensions=extensions)
        )
        shape_crs = crs_sniffer(vector_file)

        try:
            projection = CRS.from_epsg(projected_crs)
            if projection.is_geographic:
                msg = f"Desired CRS {projection.to_epsg()} is geographic. Areal analysis values will be in decimal-degree units."
                LOGGER.warning(msg)
        except Exception as e:
            msg = f"{e}: Failed to parse CRS definition. Exiting."
            LOGGER.error(msg)
            raise Exception(msg)

        # TODO: It would be good to one day refactor this to make use of RavenPy utils and gis utilities
        properties = []
        try:
            for i, layer_name in enumerate(fiona.listlayers(vector_file)):
                with fiona.open(vector_file, "r", crs=shape_crs, layer=i) as src:
                    for feature in src:
                        geom = shape(feature["geometry"])

                        multipolygon_check(geom)

                        transformed = geom_transform(
                            geom, source_crs=shape_crs, target_crs=projection
                        )
                        prop = {"id": feature["id"]}
                        prop.update(feature["properties"])
                        prop.update(geom_prop(transformed))

                        # Recompute the centroid location using the original projection
                        prop["centroid"] = geom_prop(geom)["centroid"]

                        properties.append(prop)

        except Exception as e:
            msg = f"{e}: Failed to extract features from shape {vector_file}."
            LOGGER.error(msg)
            raise Exception(msg)

        response.outputs["properties"].data = json.dumps(properties)

        return response
