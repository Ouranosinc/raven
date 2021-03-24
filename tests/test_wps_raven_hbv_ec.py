import datetime as dt

import numpy as np
import pytest
import xarray as xr
from pywps import Service
from pywps.tests import assert_response_success
from ravenpy.utilities.testdata import get_local_testdata

from raven.processes import RavenHBVECProcess

from .common import CFG_FILE, client_for, get_output, urlretrieve


class TestRavenHBVECProcess:
    def test_simple(self):
        client = client_for(
            Service(
                processes=[
                    RavenHBVECProcess(),
                ],
                cfgfiles=CFG_FILE,
            )
        )

        params = (
            "0.05984519, 4.072232, 2.001574, 0.03473693, 0.09985144, 0.5060520, 3.438486, 38.32455, "
            "0.4606565, 0.06303738, 2.277781, 4.873686, 0.5718813, 0.04505643, 0.877607, 18.94145,  "
            "2.036937, 0.4452843, 0.6771759, 1.141608, 1.024278"
        )

        datainputs = (
            "ts=files@xlink:href=file://{ts};"
            "params={params};"
            "start_date={start_date};"
            "end_date={end_date};"
            "init={init};"
            "name={name};"
            "run_name={run_name};"
            "area={area};"
            "latitude={latitude};"
            "longitude={longitude};"
            "elevation={elevation};".format(
                ts=get_local_testdata(
                    "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc",
                ),
                params=params,
                start_date=dt.datetime(2000, 1, 1),
                end_date=dt.datetime(2002, 1, 1),
                init="155,455",
                name="Salmon",
                run_name="test-hbv-ec",
                area="4250.6",
                elevation="843.0",
                latitude=54.4848,
                longitude=-123.3659,
            )
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="raven-hbv-ec",
            datainputs=datainputs,
        )

        assert_response_success(resp)
        out = get_output(resp.xml)
        assert "diagnostics" in out
        tmp_file, _ = urlretrieve(out["diagnostics"])
        tmp_content = open(tmp_file).readlines()

        # checking correctness of NSE (full period 1954-2011 would be NSE=0.591707 as template in Wiki)
        assert "DIAG_NASH_SUTCLIFFE" in tmp_content[0]
        idx_diag = tmp_content[0].split(",").index("DIAG_NASH_SUTCLIFFE")
        diag = float(tmp_content[1].split(",")[idx_diag])
        np.testing.assert_almost_equal(
            diag, 0.0186633, 4, err_msg="NSE is not matching expected value"
        )

        # checking correctness of RMSE (full period 1954-2011 would be RMSE=30.0535 as template in Wiki)
        assert "DIAG_RMSE" in tmp_content[0]
        idx_diag = tmp_content[0].split(",").index("DIAG_RMSE")
        diag = float(tmp_content[1].split(",")[idx_diag])
        np.testing.assert_almost_equal(
            diag, 35.5654, 3, err_msg="RMSE is not matching expected value"
        )

    def test_parallel(self):
        client = client_for(
            Service(
                processes=[
                    RavenHBVECProcess(),
                ],
                cfgfiles=CFG_FILE,
            )
        )

        params1 = (
            "0.05984519, 4.072232, 2.001574, 0.03473693, 0.09985144, 0.5060520, 3.438486, 38.32455, "
            "0.4606565, 0.06303738, 2.277781, 4.873686, 0.5718813, 0.04505643, 0.877607, 18.94145,  "
            "2.036937, 0.4452843, 0.6771759, 1.141608, 1.024278"
        )
        params2 = (
            "0.05, 4.07, 2.00, 0.03, 0.099, 0.506, 3.43, 38., "
            "0.47, 0.06, 2.2, 4.87, 0.5, 0.0, 0.8, 18.5,  "
            "2.0, 0.445, 0.677, 1.14, 1.02"
        )
        params3 = (
            "0.05984519, 4.072232, 2.001574, 0.03473693, 0.09985144, 0.5060520, 3.438486, 38.32455, "
            "0.4606565, 0.06303738, 2.277781, 4.873686, 0.5718813, 0.04505643, 0.877607, 18.94145,  "
            "2.036937, 0.4452843, 0.6771759, 1.141608, 1.024278"
        )

        # NOTE THAT ALL PARAMETERS IN SETS 1 and 3 ARE THE SAME. THEY SHOULD RETURN THE SAME RESULTS.

        datainputs = (
            "ts=files@xlink:href=file://{ts};"
            "params={params1};"
            "params={params2};"
            "params={params3};"
            "start_date={start_date};"
            "end_date={end_date};"
            "name={name};"
            "run_name={run_name};"
            "area={area};"
            "latitude={latitude};"
            "longitude={longitude};"
            "elevation={elevation};".format(
                ts=get_local_testdata(
                    "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc",
                ),
                params1=params1,
                params2=params2,
                params3=params3,
                start_date=dt.datetime(2000, 1, 1),
                end_date=dt.datetime(2002, 1, 1),
                name="Salmon",
                run_name="test",
                area="4250.6",
                elevation="843.0",
                latitude=54.4848,
                longitude=-123.3659,
            )
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="raven-hbv-ec",
            datainputs=datainputs,
        )

        assert_response_success(resp)
        tmp_file, _ = urlretrieve(get_output(resp.xml)["hydrograph"])
        ds = xr.open_dataset(tmp_file)

        assert ds.variables["q_sim"].shape[0] == 3

        # THIS TEST FAILS. PARAMETERS ARE NOT PASSED IN ORDER...?
        assert all(ds.variables["q_sim"][0, :] == ds.variables["q_sim"][2, :])
