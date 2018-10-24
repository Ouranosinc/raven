from raven.gis import RavenShape
import os
import glob


class TestGISClasses:
    pass

# data_dir = os.path.join(os.getcwd(), 'example_data')
# files = glob.glob(os.path.join(data_dir, '*'), recursive=True)
# shapes = [shp for shp in files if shp.endswith('.gml') or shp.endswith('.shp')]
# dems = [dem for dem in files if dem.endswith('.tif') or dem.endswith('.tiff')]
#
# print(shapes)
# print(dems)
#
# rvn = RavenShape(shapes[0])
# clip_gen = rvn.clip(dems[0])
#
# for clip in clip_gen:
#     print(type(clip))
#
# averaged_parcels = rvn.average(dems[0])
#
# for region in averaged_parcels:
#     print(region)
#
# for n in shapes:
#     rvn = RavenShape(n)
#     rvn.stats()
#     print(rvn.feature_centroids())
#     print('\n')
