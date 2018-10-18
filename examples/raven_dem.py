from osgeo import gdal
from osgeo import ogr
# import geopandas as gpd
import os
import glob


class RavenDEM(object):

    def __init__(self, dem, shape):
        self._dem = gdal.Open(dem)
        self._shape = ogr.Open(shape)
        self.__shape_filename__ = str(os.path.basename(shape))
        self.__dem_filename__ = str(os.path.basename(dem))

    def shape_stats(self):
        print('Shape stats:',  self.__shape_filename__)
        layer = self._shape.GetLayer()
        layer_def = layer.GetLayerDefn()
        print('Layer_name: ', layer_def.GetName())

        features = layer.GetFeatureCount()
        print('Number of features:', features)
        print('Field names:')
        for i in range(layer_def.GetFieldCount()):
            print(layer_def.GetFieldDefn(i).GetName())

        feature = layer.GetNextFeature()
        geo_ref = feature.GetGeometryRef()
        ref = geo_ref.GetSpatialReference()
        print('Projection: EPSG', ref.GetAuthorityCode('GEOGCS'), '\n')

    def dem_stats(self):
        print('DEM stats: ', self.__dem_filename__)

    def clip_raster_with_shape(self):
        # band = 1
        # grid = self._dem
        # shape = self._shape
        pass

    def get_feature_centroids(self):
        layer = self._shape.GetLayer()
        centroids = []
        for feature in layer:
            geom = feature.GetGeometryRef()
            centroids.append(geom.Centroid().ExportToWkt())
        return centroids

    def elevation_average(self, band=1):
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
        # rvn.dem_stats()
        print(rvn.get_feature_centroids(), '\n')

