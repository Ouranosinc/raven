from osgeo import gdal, ogr, osr
# import geopandas as gpd
import os
import glob

gdal.UseExceptions()


class RavenShape:

    @staticmethod
    def _get_crs(shape):
        """Boilerplate for quickly parsing the CRS of given shape layer."""
        ogr_shape = ogr.Open(shape)
        layer = ogr_shape.GetLayer()
        spa_ref = layer.GetSpatialRef()
        wkt = spa_ref.ExportToWkt()

        srs = osr.SpatialReference()
        srs.ImportFromWkt(wkt)

        return srs.ExportToPrettyWkt()

    def __init__(self, shape):
        self._shape = ogr.Open(shape)
        self._shape_filename = os.path.basename(shape)
        self._crs = self._get_crs(shape)

    def _get_fields(self):
        """Returns the field names from attributes table of shape layer."""
        layer = self._shape.GetLayer()
        layer_def = layer.GetLayerDefn()
        fields = []
        for i in range(layer_def.GetFieldCount()):
            fields.append(layer_def.GetFieldDefn(i).GetName())
        return fields

    def _get_layer_name(self):
        """Returns the name of the shape layer. Can lead to surprising results with ESRI Shapefiles."""
        layer = self._shape.GetLayer()
        layer_def = layer.GetLayerDefn()
        return layer_def.GetName()

    def _get_layer_bounds(self):
        """Returns the lat,lon boundaries for union of all features in layer"""
        layer = self._shape.GetLayer()
        return layer.GetExtent()

    def _get_feature_count(self):
        """Returns count of all features in layer"""
        layer = self._shape.GetLayer()
        return layer.GetFeatureCount()

    def _get_feature_bounds(self):
        raise NotImplementedError

    def _check_crs(self, crs=None, dem=None):
        """
        Ensures that the CRS of the shapefile and the specified CRS from EPSG string or DEM are compatible.
        Needed to circumvent problems that can arise from OGC-compliant WKT and ESRI-WKT CRS strings.
        Not sure of a better way to handle this issue.
        """

        if (dem and crs is None) or (crs and dem is None):
            try:
                field_srs = osr.SpatialReference()
                rvn_srs = osr.SpatialReference()

                if crs is None:
                    field = gdal.Open(dem, gdal.GA_ReadOnly)
                    field_srs.ImportFromWkt(field.GetProjection())
                else:
                    field_srs.ImportFromEPSG(crs)

                rvn_srs.ImportFromWkt(self._crs)

                if field_srs.IsGeographic() and rvn_srs.IsGeographic():
                    projcs = False
                elif not field_srs.IsGeographic() or not rvn_srs.IsGeographic():
                    projcs = True
                else:
                    # TODO: Create a check that cascades to a reprojection for clipping
                    raise Exception('Mixed geographic projection types.')

                if projcs is True:
                    structure = 'PROJCS|GEOGCS|AUTHORITY'
                else:
                    structure = 'GEOGCS|AUTHORITY'

                dem_code = field_srs.GetAttrValue(structure, 1)
                rvn_code = rvn_srs.GetAttrValue(structure, 1)

                if dem_code == rvn_code:
                    return True
                else:
                    return False

            except Exception as e:
                msg = 'Failed to read CRS/DEM: {}'.format(e)
                raise Exception(msg)

        elif crs is not None and dem is not None:
            raise KeyError('Please specify a CRS or a DEM, not both.')

        else:
            raise KeyError('Please specify a CRS or a DEM')

    def _dem_check(self, dem):
        """Call to CRS checker for entered DEMs"""
        if dem is None:
            raise Exception('No DEM specified.')
        elif isinstance(dem, str):
            self._check_crs(dem=dem)

    def stats(self):
        """Outputs stats based on shapefile contents."""
        msg = 'Shape stats: {}'.format(self._shape_filename)
        print(msg)
        print('~' * len(msg))
        print('Layer name: {}'.format(self._get_layer_name()))
        print('Layer bounds: {}'.format(self._get_layer_bounds()))
        print('Number of features: {}'.format(self._get_feature_count()))
        print('Field names: {}'.format(self._get_fields()))
        print('CRS:\n\t{}'.format(self._crs))

    def feature_centroids(self):
        """Returns the centroids respecting all features of the shape layer."""
        layer = self._shape.GetLayer()
        centroids = []
        for feature in layer:
            geom = feature.GetGeometryRef()
            centroids.append(geom.Centroid().ExportToWkt())
        return centroids

    def clip(self, dem=None):

        self._dem_check(dem)

        raster = gdal.Open(dem)
        transform = raster.GetGeoTransform()
        band = raster.GetRasterBand(1)
        array = band.ReadAsArray()

        print('Geographic Transformation:{}'.format(transform), '\n',
              'Min value:{}'.format(array.min()), '\n',
              'Mean value:{}'.format(int(array.mean())), '\n',
              'Max value:{}'.format(array.max()))

        # Courtesy of https://gis.stackexchange.com/a/195208/65343
        # shape = self._shape
        # layer = shape.GetLayer()
        #
        # if layer is not None:
        #     for feat in layer:
        #         geom = feat.GetGeometryRef()
        #         # print(geom)
        #         mx, my = geom.GetX(), geom.GetY()
        #         # print(mx, my)
        #
        #         px = int((mx - transform[0]) / transform[1])
        #         py = int((mx - transform[3]) / transform[5])
        #
        #         intval = band.ReadAsArray(px, py)#, 1, 1)
        #         print(intval)
        # else:
        #     raise Exception('Layer unable to be loaded.')
        #
        # return intval
        raise NotImplementedError

    def average_elev(self, dem=None):
        pass


if __name__ == '__main__':
    gdal.GetDriverByName('GTiff')

    data_dir = os.path.join(os.getcwd(), 'example_data')
    files = glob.glob(os.path.join(data_dir, '*'), recursive=True)
    shapes = [shp for shp in files if shp.endswith('.gml') or shp.endswith('.shp')]
    dems = [dem for dem in files if dem.endswith('.tif') or dem.endswith('.tiff')]

    print(shapes)
    print(dems)

    rvn = RavenShape(shapes[0])
    output = rvn.clip(dems[0])

    print(output)

    # for n in shapes:
    #     rvn = RavenShape(n)
    #     rvn.stats()
    #     print(rvn._check_crs(crs=4326))
    #     print(rvn._check_crs(dem=dem))
    #     print(rvn.feature_centroids())
    #     print('\n')
