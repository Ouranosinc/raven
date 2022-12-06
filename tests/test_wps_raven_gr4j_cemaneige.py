import datetime as dt
import json
import zipfile
from dataclasses import asdict
from urllib.request import urlretrieve

import numpy as np
import xarray as xr
from pywps import Service
from pywps.tests import assert_response_success
from ravenpy.models import GR4JCN

from raven.processes import RavenGR4JCemaNeigeProcess

from .common import CFG_FILE, client_for, get_output


class TestRavenGR4JCemaNeigeProcess:
    def test_simple(self, get_local_testdata):
        client = client_for(
            Service(processes=[RavenGR4JCemaNeigeProcess()], cfgfiles=CFG_FILE)
        )

        params = "0.529, -3.396, 407.29, 1.072, 16.9, 0.947"

        # some params in Raven input files are derived from those 21 parameters
        # pdefaults.update({'GR4J_X1_hlf':            pdefaults['GR4J_X1']*1000./2.0})    --> x1 * 1000. / 2.0
        # pdefaults.update({'one_minus_CEMANEIGE_X2': 1.0 - pdefaults['CEMANEIGE_X2']})   --> 1.0 - x6

        datainputs = (
            "ts=files@xlink:href=file://{ts};"
            "params={params};"
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
                params=params,
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
            identifier="raven-gr4j-cemaneige",
            datainputs=datainputs,
        )

        assert_response_success(resp)

        out = get_output(resp.xml)
        assert "diagnostics" in out
        tmp_file, _ = urlretrieve(out["diagnostics"])
        tmp_content = open(tmp_file).readlines()

        # checking correctness of NSE (full period 1954-2010 would be NSE=0.511214)
        assert "DIAG_NASH_SUTCLIFFE" in tmp_content[0]
        idx_diag = tmp_content[0].split(",").index("DIAG_NASH_SUTCLIFFE")
        diag = float(tmp_content[1].split(",")[idx_diag])
        np.testing.assert_almost_equal(
            diag, -0.117301, 4, err_msg="NSE is not matching expected value"
        )

        # checking correctness of RMSE (full period 1954-2010 would be RMSE=32.8827)
        assert "DIAG_RMSE" in tmp_content[0]
        idx_diag = tmp_content[0].split(",").index("DIAG_RMSE")
        diag = float(tmp_content[1].split(",")[idx_diag])
        np.testing.assert_almost_equal(
            diag, 37.9493, 4, err_msg="RMSE is not matching expected value"
        )

        assert "rv_config" in out
        rv_config, _ = urlretrieve(out["rv_config"])
        z = zipfile.ZipFile(rv_config)
        assert len(z.filelist) == 5

    def test_simple_alternate_hrus_param(self, get_local_testdata):
        client = client_for(
            Service(processes=[RavenGR4JCemaNeigeProcess()], cfgfiles=CFG_FILE)
        )

        params = "0.529, -3.396, 407.29, 1.072, 16.9, 0.947"

        # some params in Raven input files are derived from those 21 parameters
        # pdefaults.update({'GR4J_X1_hlf':            pdefaults['GR4J_X1']*1000./2.0})    --> x1 * 1000. / 2.0
        # pdefaults.update({'one_minus_CEMANEIGE_X2': 1.0 - pdefaults['CEMANEIGE_X2']})   --> 1.0 - x6

        salmon_land_hru_1 = dict(
            area=4250.6, elevation=843.0, latitude=54.4848, longitude=-123.3659
        )
        hrus = (GR4JCN.LandHRU(**salmon_land_hru_1),)
        hrus = json.dumps(list(map(asdict, hrus)))

        datainputs = (
            "ts=files@xlink:href=file://{ts};"
            "params={params};"
            "start_date={start_date};"
            "end_date={end_date};"
            "name={name};"
            "run_name={run_name};"
            "hrus={hrus};".format(
                ts=get_local_testdata(
                    "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc",
                ),
                params=params,
                start_date=dt.datetime(2000, 1, 1),
                end_date=dt.datetime(2002, 1, 1),
                name="Salmon",
                run_name="test",
                hrus=hrus,
            )
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="raven-gr4j-cemaneige",
            datainputs=datainputs,
        )

        assert_response_success(resp)

        out = get_output(resp.xml)
        assert "diagnostics" in out
        tmp_file, _ = urlretrieve(out["diagnostics"])
        tmp_content = open(tmp_file).readlines()

        # checking correctness of NSE (full period 1954-2010 would be NSE=0.511214)
        assert "DIAG_NASH_SUTCLIFFE" in tmp_content[0]
        idx_diag = tmp_content[0].split(",").index("DIAG_NASH_SUTCLIFFE")
        diag = float(tmp_content[1].split(",")[idx_diag])
        np.testing.assert_almost_equal(
            diag, -0.117301, 4, err_msg="NSE is not matching expected value"
        )

        # checking correctness of RMSE (full period 1954-2010 would be RMSE=32.8827)
        assert "DIAG_RMSE" in tmp_content[0]
        idx_diag = tmp_content[0].split(",").index("DIAG_RMSE")
        diag = float(tmp_content[1].split(",")[idx_diag])
        np.testing.assert_almost_equal(
            diag, 37.9493, 4, err_msg="RMSE is not matching expected value"
        )

        assert "rv_config" in out
        rv_config, _ = urlretrieve(out["rv_config"])
        z = zipfile.ZipFile(rv_config)
        assert len(z.filelist) == 5

    def test_parallel(self, get_local_testdata):
        client = client_for(
            Service(processes=[RavenGR4JCemaNeigeProcess()], cfgfiles=CFG_FILE)
        )

        params1 = "0.529, -3.396, 407.29, 1.072, 16.9, 0.947"
        params2 = "0.5, -3.3, 407.2, 1.0, 16.6, 0.9"

        # some params in Raven input files are derived from those 21 parameters
        # pdefaults.update({'GR4J_X1_hlf':            pdefaults['GR4J_X1']*1000./2.0})    --> x1 * 1000. / 2.0
        # pdefaults.update({'one_minus_CEMANEIGE_X2': 1.0 - pdefaults['CEMANEIGE_X2']})   --> 1.0 - x6

        datainputs = (
            "ts=files@xlink:href=file://{ts};"
            "params={params1};"
            "params={params2};"
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
            identifier="raven-gr4j-cemaneige",
            datainputs=datainputs,
        )

        assert_response_success(resp)
        tmp_file, _ = urlretrieve(get_output(resp.xml)["hydrograph"])
        ds = xr.open_dataset(tmp_file)
        assert ds.variables["q_sim"].shape[0] == 2

    def test_ravenpy_error(self, get_local_testdata):
        client = client_for(
            Service(processes=[RavenGR4JCemaNeigeProcess()], cfgfiles=CFG_FILE)
        )

        # Deliberate error: no HRUs!

        datainputs = (
            "ts=files@xlink:href=file://{ts};"
            "params={params};"
            "start_date={start_date};"
            "end_date={end_date};"
            "name={name};"
            "run_name={run_name};".format(
                ts=get_local_testdata(
                    "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc",
                ),
                params="0.529, -3.396, 407.29, 1.072, 16.9, 0.947",
                start_date=dt.datetime(2000, 1, 1),
                end_date=dt.datetime(2002, 1, 1),
                name="Salmon",
                run_name="test",
            )
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="raven-gr4j-cemaneige",
            datainputs=datainputs,
        )

        assert "CHydroUnit constructor:: HRU 1 has a negative or zero area" in str(
            resp.data
        )
