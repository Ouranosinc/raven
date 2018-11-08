from raven.gis import RavenShape, RavenRaster
import os
# import pytest


class TestGISClasses:

    def test_raven_shape(self):
        here = os.getcwd()

        dem = os.path.join(here, 'tests/testdata/earthenv_dem_90m/earthenv_dem90_southernQuebec.tiff')
        shape = os.path.join(here, 'tests/testdata/donneesqc_mrc_poly/donnees_quebec_mrc_polygones.gml')

        rvn = RavenShape(shape)
        assert rvn.crs() == '4326'
        assert not rvn.check_crs(crs=666)
        assert rvn.check_crs(raster=dem)
        assert len(rvn) == 137
        assert rvn.clip(dem)

    def test_raven_raster(self):
        here = os.getcwd()

        dem = os.path.join(here, 'tests/testdata/earthenv_dem_90m/earthenv_dem90_southernQuebec.tiff')
        shape = os.path.join(here, 'tests/testdata/donneesqc_mrc_poly/donnees_quebec_mrc_polygones.gml')

        rvn = RavenRaster(dem)
        assert rvn.crs() == '4326'
        assert not rvn.check_crs(crs='5555')
        assert rvn.check_crs(shape=shape)


# if __name__ == "__main__":
#     here = os.getcwd()
#
#     dem = os.path.join(here, '../tests/testdata/earthenv_dem_90m/earthenv_dem90_southernQuebec.tiff')
#     shape = os.path.join(here, '../tests/testdata/donneesqc_mrc_poly/donnees_quebec_mrc_polygones.gml')
#
#     rvn = RavenShape(shape)
#
#     print(rvn)
#
#     print(rvn.layer_bounds())
#
#     for parcel in rvn.feature_bounds():
#         print(parcel)
#
#     for parcel in rvn.area(crs=3857):
#         print(parcel)
#
#     for centroid in rvn.centroids():
#         print(centroid)
#
#     for clip in rvn.clip(dem):
#         print(clip)
#
#     for region in rvn.average(raster=dem):
#         print(region)
#
#     # rvn = RavenRaster(dem)
#
#     rvn.warp(crs=4269)
#     print(rvn.projected())
#
#     print(rvn.pixel_counts())
