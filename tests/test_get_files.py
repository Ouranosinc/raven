from raven.tutorial import get_file, open_dataset, _default_cache_dir
import xarray
from pathlib import Path

git_url = "https://github.com/Ouranosinc/raven-testdata"


class TestRemoteFileAccess:
    git_url = "https://github.com/Ouranosinc/raven-testdata"
    branch = "master"
    dap_url = 'http://test.opendap.org:80/opendap/data/nc/'

    def test_get_file(self):
        file = get_file(name="ostrich-hbv-ec/raven-hbv-salmon.rvi", branch=self.branch)
        assert Path(_default_cache_dir).exists()
        assert file.is_file()
        with file.open() as f:
            header = f.read()
            assert ":FileType          rvi ASCII Raven 2.8.2" in header

    def test_open_dataset(self):
        ds = open_dataset(
            name="raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc",
            branch=self.branch,
        )
        assert (
            Path(_default_cache_dir)
            .joinpath(
                "raven-gr4j-cemaneige", "Salmon-River-Near-Prince-George_meteo_daily.nc"
            )
            .exists()
        )
        assert isinstance(ds, xarray.Dataset)

    def test_open_dataset_no_cache(self):
        ds = open_dataset(
            name="raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily_3d.nc",
            branch=self.branch,
            cache=False
        )
        assert (
            not Path(_default_cache_dir)
            .joinpath(
                "raven-gr4j-cemaneige",
                "Salmon-River-Near-Prince-George_meteo_daily_3d.nc",
            )
            .exists()
        )
        assert isinstance(ds, xarray.Dataset)

    def test_dap_access(self):
        ds = open_dataset(
            name="20070917-MODIS_A-JPL-L2P-A2007260000000.L2_LAC_GHRSST-v01.nc",
            dap_url=self.dap_url
        )
        assert isinstance(ds, xarray.Dataset)
