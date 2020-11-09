import xarray as xr
import xclim.sdba as sdba
from xclim import subset


class TestBiasCorrect:
    def test_bias_correction(self):

        fut_data = (
            "https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/datasets/simulations/bias_adjusted/"
            "cmip5/nasa/nex-gddp-1.0/day_inmcm4_historical+rcp85_nex-gddp.ncml"
        )
        ref_data = (
            "https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/datasets/simulations/bias_adjusted/"
            "cmip5/nasa/nex-gddp-1.0/day_inmcm4_historical+rcp45_nex-gddp.ncml"
        )
        hist_data = "https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/datasets/gridded_obs/nrcan_v2.ncml"

        lat = 49.5
        lon = -72.96

        # Open these datasets
        ds_fut = xr.open_dataset(fut_data)
        ds_ref = xr.open_dataset(ref_data)
        ds_his = xr.open_dataset(hist_data)

        # Subset the data to the desired location (2x2 degree box)
        ds_fut_sub = subset.subset_bbox(
            ds_fut,
            lon_bnds=[lon - 1, lon + 1],
            lat_bnds=[lat - 1, lat + 1],
            start_date="2070",
            end_date="2071",
        ).mean(dim={"lat", "lon"}, keep_attrs=True)
        ds_ref_sub = subset.subset_bbox(
            ds_ref,
            lon_bnds=[lon - 1, lon + 1],
            lat_bnds=[lat - 1, lat + 1],
            start_date="1971",
            end_date="1972",
        ).mean(dim={"lat", "lon"}, keep_attrs=True)
        ds_his_sub = subset.subset_bbox(
            ds_his,
            lon_bnds=[lon - 1, lon + 1],
            lat_bnds=[lat - 1, lat + 1],
            start_date="1971",
            end_date="1972",
        ).mean(dim={"lat", "lon"}, keep_attrs=True)

        group_month_nowindow = sdba.utils.Grouper("time.month")
        Adj = sdba.DetrendedQuantileMapping(
            nquantiles=50, kind="+", group=group_month_nowindow
        )
        # Train the model to find the correction factors
        Adj.train(ds_ref_sub["pr"], ds_his_sub["pr"])

        # Apply the factors to the future data to bias-correct
        Adj.adjust(ds_fut_sub["pr"], interp="linear")

        # Repeat for temperature max
        Adj.train(ds_ref_sub["tasmax"], ds_his_sub["tasmax"])

        # Apply the factors to the future data to bias-correct
        Adj.adjust(ds_fut_sub["tasmax"], interp="linear")

        # Repeat for tasmin
        Adj.train(ds_ref_sub["tasmin"], ds_his_sub["tasmin"])
        Adj.adjust(ds_fut_sub["tasmin"], interp="linear")

        # TODO: Add numerical check
