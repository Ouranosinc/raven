from raven.gis import RavenShape
import os
import glob


class TestGISClasses:

    def test_raven_shape(self):
        pass

    def test_raven_raster(self):
        pass

    def test_raven_gis(self):
        pass

# shapes = "/home/tjs/git/raven/tests/testdata/donneesqc_mrc_poly/donnees_quebec_mrc_polygones.gml"
#
# dems = "/home/tjs/git/raven/tests/testdata/earthenv_dem_90m/earthenv_dem90_southernQuebec.tiff"
#
# print(shapes)
# print(dems)
#
# rvn = RavenShape(shapes)
# rvn.stats()
# clip_gen = rvn.clip(dems)
#
# for clip in clip_gen:
#     print(type(clip))
#     print(clip)
#
# averaged_parcels = rvn.average(dems)
#
# for region in averaged_parcels:
#     print(region)
#
#
# rvn = RavenShape(shapes)
# rvn.stats()
# print(rvn.feature_centroids())
# print('\n')
