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


# print(shapes)
# print(dems)
#
# rvn = RavenShape(shapes)

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
# print(rvn.feature_centroids())
# print('\n')
