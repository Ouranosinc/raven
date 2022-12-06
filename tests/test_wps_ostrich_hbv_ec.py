import datetime as dt
from urllib.request import urlretrieve

import numpy as np
from pywps import Service
from pywps.tests import assert_response_success

from raven.processes import OstrichHBVECProcess

from .common import CFG_FILE, client_for, get_output


class TestOstrichHBVECProcess:
    def test_simple(self, get_local_testdata):
        client = client_for(
            Service(processes=[OstrichHBVECProcess()], cfgfiles=CFG_FILE)
        )

        lowerBounds = "-3.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.01, 0.05, 0.01, \
                      0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.05, 0.8, 0.8"
        upperBounds = "3.0, 8.0, 8.0, 0.1, 1.0, 1.0, 7.0, 100.0, 1.0, 0.1, 6.0, 5.0, \
                      5.0, 0.2, 1.0, 30.0, 3.0, 2.0, 1.0, 1.5, 1.5"

        datainputs = (
            "ts=files@xlink:href=file://{ts};"
            "algorithm={algorithm};"
            "max_iterations={max_iterations};"
            "lowerBounds={lowerBounds};"
            "upperBounds={upperBounds};"
            "start_date={start_date};"
            "duration={duration};"
            "name={name};"
            "run_name={run_name};"
            "area={area};"
            "latitude={latitude};"
            "longitude={longitude};"
            "elevation={elevation};"
            "random_numbers=file@xlink:href=file://{random_numbers};"
            "random_seed=0".format(
                ts=get_local_testdata(
                    "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc",
                ),
                algorithm="DDS",
                max_iterations=10,
                lowerBounds=lowerBounds,
                upperBounds=upperBounds,
                start_date=dt.datetime(1954, 1, 1),
                duration=208,
                name="Salmon",
                run_name="test",
                area="4250.6",
                elevation="843.0",
                latitude=54.4848,
                longitude=-123.3659,
                random_numbers=get_local_testdata(
                    "ostrich-gr4j-cemaneige/OstRandomNumbers.txt"
                ),
            )
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="ostrich-hbv-ec",
            datainputs=datainputs,
        )

        assert_response_success(resp)

        out = get_output(resp.xml)
        assert "diagnostics" in out
        tmp_file, _ = urlretrieve(out["diagnostics"])
        tmp_content = open(tmp_file).readlines()

        # TODO Julie :: values not adjusted yet!!! WPS needs to work first ...

        # checking correctness of NSE (full period 1954-2010 with budget 50 would be NSE=6.034670E-01)
        assert "DIAG_NASH_SUTCLIFFE" in tmp_content[0]
        idx_diag = tmp_content[0].split(",").index("DIAG_NASH_SUTCLIFFE")
        diag = float(tmp_content[1].split(",")[idx_diag])
        np.testing.assert_almost_equal(
            diag, -0.225991, 4, err_msg="NSE is not matching expected value"
        )

        # checking correctness of RMSE (full period 1954-2010 would be RMSE=????)
        assert "DIAG_RMSE" in tmp_content[0]
        idx_diag = tmp_content[0].split(",").index("DIAG_RMSE")
        diag = float(tmp_content[1].split(",")[idx_diag])
        np.testing.assert_almost_equal(
            diag, 57.3687, 4, err_msg="RMSE is not matching expected value"
        )
