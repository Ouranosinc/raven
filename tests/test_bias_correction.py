import datetime as dt

import xarray as xr
import xclim.sdba as sdba
from ravenpy.utilities.testdata import get_local_testdata
from xclim import subset


class TestBiasCorrect:
    def test_bias_correction(self):

        ds_fut_sub = xr.open_dataset(
            get_local_testdata(
                "cmip5/nasa_nex-gddp-1.0_day_inmcm4_historical+rcp85_nex-gddp_2070-2071_subset.nc",
            )
        )
        ds_ref_sub = xr.open_dataset(
            get_local_testdata(
                "cmip5/nasa_nex-gddp-1.0_day_inmcm4_historical+rcp45_nex-gddp_1971-1972_subset.nc",
            )
        )
        ds_his_sub = xr.open_dataset(
            get_local_testdata("nrcan/NRCAN_1971-1972_subset.nc")
        )

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
