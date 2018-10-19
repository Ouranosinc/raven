from osgeo import gdal
from osgeo import ogr
# import geopandas as gpd
import os
import glob


class RavenDEM(object):

    def _get_proj(self):
        layer = self._shape.GetLayer()
        feature = layer.GetNextFeature()
        geo_ref = feature.GetGeometryRef()
        ref = geo_ref.GetSpatialReference()
        proj = 'EPSG:{}'.format(ref.GetAuthorityCode("GEOGCS"))
        return {'Projection': proj}

    def _get_fields(self):
        layer = self._shape.GetLayer()
        layer_def = layer.GetLayerDefn()
        fields = []
        for i in range(layer_def.GetFieldCount()):
            fields.append(layer_def.GetFieldDefn(i).GetName())
        return fields

    def _get_layer(self):
        layer = self._shape.GetLayer()
        layer_def = layer.GetLayerDefn()
        return layer_def.GetName()

    def _get_features(self):
        layer = self._shape.GetLayer()
        return layer.GetFeatureCount()

    def __init__(self, dem, shape):
        self._dem = gdal.Open(dem)
        self._shape = ogr.Open(shape)
        self.__shape_filename__ = str(os.path.basename(shape))
        self.__dem_filename__ = str(os.path.basename(dem))
        self.__crs__ = self._get_proj()

    def shape_stats(self):
        print('Shape stats:',  self.__shape_filename__)
        print('Layer_name: ', self._get_layer())
        print('Number of features:', self._get_features())
        print('Field names:', self._get_fields())
        print('CRS:', self._get_proj())

    def dem_stats(self):
        print('DEM stats: ', self.__dem_filename__)

    def reproject(self, crs):
        pass

    def clip_raster_with_shape(self):
        # band = 1
        # grid = self._dem
        # shape = self._shape
        pass

    def feature_centroids(self):
        layer = self._shape.GetLayer()
        centroids = []
        for feature in layer:
            geom = feature.GetGeometryRef()
            centroids.append(geom.Centroid().ExportToWkt())
        return centroids

    def average_elev(self, band=1):
        # band = 1
        pass


if __name__ == '__main__':

    data_dir = os.path.join(os.getcwd(), 'example_data')
    files = glob.glob(os.path.join(data_dir, '*'), recursive=True)
    shapes = [shp for shp in files if shp.endswith('.gml') or shp.endswith('.shp')]
    dems = [dem for dem in files if dem.endswith('.tif') or dem.endswith('.tiff')]

    print(shapes)

    for n in shapes:
        rvn = RavenDEM(dems[0], n)
        rvn.shape_stats()
        print(rvn.feature_centroids(), '\n')

