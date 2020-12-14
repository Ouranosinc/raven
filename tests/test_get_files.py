from pathlib import Path

import xarray

from raven.tutorial import _default_cache_dir, get_file, open_dataset




class TestRemoteFileAccess:
    dap_url = "http://test.opendap.org:80/opendap/data/nc/"
    git_url = "https://github.com/Ouranosinc/raven-testdata"
    branch = "master"
    DCD = Path(_default_cache_dir) / branch

    def test_get_file(self):
        file = get_file(name="ostrich-hbv-ec/raven-hbv-salmon.rvi", branch=self.branch)

        assert self.DCD.exists()
        assert file.is_file()
        with file.open() as f:
            header = f.read()
            assert ":FileType          rvi ASCII Raven 2.8.2" in header

    def test_open_dataset(self):
        ds = open_dataset(
            name="raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc",
            branch=self.branch,
        )

        assert isinstance(ds, xarray.Dataset)
        assert (self.DCD / "raven-gr4j-cemaneige" / "Salmon-River-Near-Prince-George_meteo_daily.nc").exists()

    def test_open_dataset_no_cache(self):
        ds = open_dataset(
            name="raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily_3d.nc",
            branch=self.branch,
            cache=False,
        )

        assert isinstance(ds, xarray.Dataset)
        assert not (self.DCD / "raven-gr4j-cemaneige" / "Salmon-River-Near-Prince-George_meteo_daily_3d.nc").exists()

    def test_dap_access(self):
        ds = open_dataset(
            name="20070917-MODIS_A-JPL-L2P-A2007260000000.L2_LAC_GHRSST-v01.nc",
            dap_url=self.dap_url,
        )

        assert isinstance(ds, xarray.Dataset)
