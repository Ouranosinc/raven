from osgeo import gdal, ogr
import pyproj
# import geopandas as gpd
import os
import glob


class RavenShape(object):

    def _get_proj(self):

        layer = self._shape.GetLayer()
        feature = layer.GetNextFeature()
        geo_ref = feature.GetGeometryRef()
        ref = geo_ref.GetSpatialReference()

        return ref.GetAuthorityCode("GEOGCS")

    def _get_pyprojCRS(self):

        proj = 'epsg:{}'.format(self._get_proj())
        return pyproj.Proj(init=proj)

    def _get_fields(self):

        layer = self._shape.GetLayer()
        layer_def = layer.GetLayerDefn()
        fields = []
        for i in range(layer_def.GetFieldCount()):
            fields.append(layer_def.GetFieldDefn(i).GetName())
        return fields

    def _get_layer_name(self):

        layer = self._shape.GetLayer()
        layer_def = layer.GetLayerDefn()
        return layer_def.GetName()

    def _get_layer_bounds(self):

        layer = self._shape.GetLayer()
        return layer.GetExtent()

    def _get_feature_count(self):

        layer = self._shape.GetLayer()
        return layer.GetFeatureCount()

    def _get_feature_bounds(self):
        pass

    def __init__(self, shape):
        self._shape = ogr.Open(shape)
        self.__shape_filename__ = str(os.path.basename(shape))
        self.__crs__ = self._get_pyprojCRS()

    def shape_stats(self):

        msg = 'Shape stats: {}'.format(self.__shape_filename__)
        print(msg)
        print('~' * len(msg))
        print('Layer name: ', self._get_layer_name())
        print('Layer bounds:', self._get_layer_bounds())
        print('Number of features:', self._get_feature_count())
        print('Field names:', self._get_fields())
        print('CRS: EPSG:{}'.format(self._get_proj()))

    def check_crs(self, crs=None, dem=None):

        if crs and dem:
            msg = 'Please specify a crs or a dem, not both.'
            raise KeyError(msg)

        elif crs:
            if crs == self.__crs__:
                print('Shape is already in this CRS')
                return None

            try:
                pass
            except Exception as e:
                msg = 'CRS cannot be parsed: {}'.format(e)
                raise ValueError(msg)

        elif dem:
            msg = 'Please specify a crs or a dem'
            raise KeyError(msg)

    def feature_centroids(self):

        layer = self._shape.GetLayer()
        centroids = []
        for feature in layer:
            geom = feature.GetGeometryRef()
            centroids.append(geom.Centroid().ExportToWkt())
        return centroids

    def clip_raster(self, dem=None):

        raster = gdal.Open(dem)
        # band = 1
        # grid = self._dem
        # shape = self._shape
        pass

    def average_raster_elev(self, dem=None):
        # band = 1
        pass


if __name__ == '__main__':

    data_dir = os.path.join(os.getcwd(), 'example_data')
    files = glob.glob(os.path.join(data_dir, '*'), recursive=True)
    shapes = [shp for shp in files if shp.endswith('.gml') or shp.endswith('.shp')]
    dems = [dem for dem in files if dem.endswith('.tif') or dem.endswith('.tiff')]

    print(shapes)

    for n in shapes:
        rvn = RavenShape(n)
        rvn.shape_stats()
        print(rvn.feature_centroids())
        print('\n')

