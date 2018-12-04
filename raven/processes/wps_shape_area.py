import logging

import re
import fiona
from shapely.geometry import shape
from pywps import LiteralInput, ComplexInput
from pywps import LiteralOutput
from pywps import Process

LOGGER = logging.getLogger("PYWPS")


class ShapeArea(Process):
    """Given a file containing vector data, provide general information and spatial characteristics"""

    def __init__(self):
        inputs = [
            LiteralInput('expand features', 'Examine all features in shapefile', data_type='boolean'),
            ComplexInput('shape', 'Vector Shape',
                         abstract='An URL pointing to either an ESRI Shapefile, GML, GeoJSON, or any other file in a'
                                  ' standard vector format. The ESRI Shapefile must be zipped and contain the .shp,'
                                  ' .shx, and .dbf. The shape CRS definition should also match the DEM CRS.',
                         min_occurs=1,
                         default=None, supported_formats=[])]
        # supported_formats=[
        #            # shp
        #            {mimeType: 'application/zip',
        #            encoding: '???',
        #            schema: None},
        #
        #            # gml
        #            {mimeType: 'text/xml',
        #            encoding:'utf-8',
        #            schema:'http://schemas.opengis.net/gml/3.2.1/gml.xsd'},
        #
        #            # json
        #            {mimeType: 'text/plain',
        #            encoding: 'iso-8859-2',
        #            schema: None
        #            },
        #
        #            # kml
        #            {mimeType: 'text/xml',
        #            encoding: 'windows-1250',
        #            schema: 'http://schemas.opengis.net/kml/2.2.0/ogckml22.xsd'}
        #            ])

        outputs = [
            LiteralOutput('area', 'Area Calculations', data_type='float', abstract='Area of shape in sq. km.'),
            LiteralOutput('centroids', 'Centroid Locations', data_type='float',
                          abstract="Geographic locations of feature centroids"),
            LiteralOutput('schemas', 'Feature schemas',
                          abstract='Geographic representations and descriptions of shape features')
        ]

        super(ShapeArea, self).__init__(
            self._shapearea_handler,
            identifier="shape_area",
            title="Shape Area",
            version="1.0",
            abstract="Return shape area in square Kilometres based on line boundaries of a polygonal vector file.",
            metadata=[],
            inputs=inputs,
            outputs=outputs,
            status_supported=True,
            store_supported=True)

    @staticmethod
    def _shapearea_handler(request, response):

        def address_append(address):
            zipped = re.search(r'(\.zip)', address)
            tarred = re.search(r'(\.tar)', address)

            try:
                if zipped:
                    return 'zip://{}'.format(address)
                elif tarred:
                    return 'tar://{}'.format(address)
                else:
                    return address
            except Exception as e:
                msg = 'Failed to prefix or parse URL: {}'.format(e)
                raise Exception(msg)

        shape_fn = request.inputs['shape'][0].file
        shape_url = address_append(shape_fn)

        schemas = []
        centroids = []
        areas = []

        features = False
        if request.inputs['expand features'][0]:
            features = True

        for i, layername in enumerate(fiona.listlayers(shape_url)):
            if features:
                with fiona.open(shape_url, 'r', layer=i) as src:
                    geom = shape(src['geometry'])

                    schemas.append(src.schema)
                    centroids.append(geom.centroid)
                    areas.append(geom.area)

                    src.close()
            else:
                with fiona.open(shape_url, 'r') as src:
                    geom = shape(src['geometry'])

                    schemas.append(src.schema)
                    centroids.append(geom.centroid)
                    areas.append(geom.area)

                    src.close()

        response.outputs['area'].data = areas
        response.outputs['centroids'].data = centroids
        response.outputs['schemas'].data = schemas

        return response
