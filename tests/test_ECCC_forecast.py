import datetime as dt
import tempfile

import pandas as pd
from fiona.io import ZipMemoryFile

from raven.models import GR4JCN
from raven.utilities import forecasting

from .common import TESTDATA

"""
Test to perform a hindcast using auto-queried ECCC data aggregated on THREDDS.
Currently only runs GEPS, eventually will run GEPS, GDPS, REPS and RDPS.
To do so will need to add the actual data from ECCC but this is a proof of concept.
"""


class TestECCCForecast:
    def test_forecasting_GEPS(self, tmpdir):

        # Prepare a RAVEN model run using historical data, GR4JCN in this case.
        # This is a dummy run to get initial states. In a real forecast situation,
        # this run would end on the day before the forecast, but process is the same.
        ts = TESTDATA["raven-gr4j-cemaneige-nc-ts"]
        model = GR4JCN(workdir=tmpdir)
        model(
            ts,
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2002, 6, 1),
            area=44250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            params=(0.529, -3.396, 407.29, 1.072, 16.9, 0.947),
        )

        # Extract the final states that will be used as the next initial states
        rvc = model.outputs["solution"]

        # Collect the most recent forecast data for location and climate model.
        # Limited to GEPS for now
        with open(TESTDATA["watershed_vector"], "rb") as f:
            with ZipMemoryFile(f.read()) as z:
                region_coll = z.open("LSJ_LL.shp")
                fcst = forecasting.get_recent_ECCC_forecast(
                    region_coll, climate_model="GEPS"
                )

        # write the forecast data to file
        ts = f"{tmpdir}/fcstfile.nc"
        fcst.to_netcdf(ts)

        # It is necessary to clean the model state because the input variables of the previous
        # model are not the same as the ones provided in the forecast model. therefore, if we
        # do not clean, the model will simply add the hindcast file to the list of available
        # data provided in the testdata above. Then the dates will not work, and the model errors.

        model = GR4JCN()

        model.rvc.parse(rvc.read_text())

        model(
            ts=(ts,),
            nc_index=range(fcst.dims.get("member")),
            duration=9,
            area=44250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            params=(0.529, -3.396, 407.29, 1.072, 16.9, 0.947),
            overwrite=True,
            pr={"time_shift": -0.25, "deaccumulate": True},
        )

        # The model now has the forecast data generated and it has 10 days of forecasts.
        assert len(model.q_sim.values) == 10
        # Also see if GEPS has 20 members produced.
        assert model.q_sim.values.shape[1] == 20

