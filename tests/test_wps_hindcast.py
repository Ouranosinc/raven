import datetime as dt
import json
import tempfile

import matplotlib.pyplot as plt
import pytest
import xarray as xr
from pywps import Service
from pywps.tests import assert_response_success
from ravenpy.utilities.testdata import get_local_testdata
import numpy as np

from raven.processes import HindcastingProcess

from .common import CFG_FILE, client_for, get_output, urlretrieve


@pytest.mark.online
class TestHindcasting:
    def test_GEPS(self):
        client = client_for(
            Service(
                processes=[
                    HindcastingProcess(),
                ],
                cfgfiles=CFG_FILE,
            )
        )
        #
        # model = 'HMETS'
        # params = '9.5019, 0.2774, 6.3942, 0.6884, 1.2875, 5.4134, 2.3641, 0.0973, 0.0464, 0.1998, 0.0222, -1.0919, ' \
        #         '2.6851, 0.3740, 1.0000, 0.4739, 0.0114, 0.0243, 0.0069, 310.7211, 916.1947'

        params = "0.529, -3.396, 407.29, 1.072, 16.9, 0.947"
        forecast_model = "GEPS"
        region_vector = get_local_testdata("watershed_vector/LSJ_LL.zip")
        pr = json.dumps(
            {
                "pr": {
                    "linear_transform": (1.0, 0.0),
                    "time_shift": -0.25,
                    "deaccumulate": True,
                }
            }
        )
        tas = json.dumps({"tas": {"linear_transform": (1.0, 0.0), "time_shift": -0.25}})
        rvc = get_local_testdata("gr4j_cemaneige/solution.rvc")
        hdate = dt.datetime(2018, 6, 1)

        # Date of the forecast that will be used to determine the members of the climatology-based ESP
        # (same day of year of all other years)
        datainputs = (
            "gr4jcn={params};"
            "latitude={latitude};"
            "longitude={longitude};"
            "name={name};"
            "area={area};"
            "duration={duration};"
            "elevation={elevation};"
            "forecast_model={forecast_model};"
            "region_vector=file@xlink:href=file://{region_vector};"
            "rain_snow_fraction={rain_snow_fraction};"
            "nc_spec={pr};"
            "nc_spec={tas};"
            "rvc=file@xlink:href=file://{rvc};"
            "hdate={hdate};".format(
                params=params,
                latitude=54.4848,
                longitude=-123.3659,
                name="Salmon",
                area="4250.6",
                duration=8,
                elevation="843.0",
                forecast_model=forecast_model,
                region_vector=region_vector,
                rain_snow_fraction="RAINSNOW_DINGMAN",
                pr=pr,
                tas=tas,
                rvc=rvc,
                hdate=hdate,
            )
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="hindcasting",
            datainputs=datainputs,
        )

        assert_response_success(resp)
        out = get_output(resp.xml)
        assert "hydrograph" in out
        forecast, _ = urlretrieve(out["hydrograph"])
        q_sim = xr.open_dataset(forecast)["q_sim"]
        np.testing.assert_almost_equal(q_sim.isel(time=-1).mean(), [12.585823219473196])
        assert "member" in q_sim.dims
        # To display the forecast
        # q_sim.plot(); plt.show()

