from raven.utilities import gis


class TestSelect_hybas_domain:

    def test_na(self):
        bbox = (-68., 50., -68., 50.)
        dom = gis.select_hybas_domain(bbox)
        assert dom == 'na'

    def test_ar(self):
        bbox = (-114.65, 61.35, -114.65, 61.35)
        dom = gis.select_hybas_domain(bbox)
        assert dom == 'ar'
