import datetime as dt
import json

import xarray as xr
from pywps import Service
from pywps.tests import assert_response_success
from ravenpy.utilities.testdata import get_test_data

from raven.processes import RavenHMETSProcess

from .common import CFG_FILE, client_for


class TestRavenNRCANProcess:
    def test_simple(self):
        client = client_for(
            Service(
                processes=[
                    RavenHMETSProcess(),
                ],
                cfgfiles=CFG_FILE,
            )
        )

        salmon = get_test_data(
            "raven-gr4j-cemaneige", "Salmon-River-Near-Prince-George_meteo_daily.nc"
        )[0]
        salmon = xr.open_dataset(salmon)
        lat = salmon.lat.values[0]
        lon = salmon.lon.values[0]

        ts = get_test_data("nrcan", "NRCAN_ts.nc")[0]

        params = (
            "9.5019, 0.2774, 6.3942, 0.6884, 1.2875, 5.4134, 2.3641, 0.0973, 0.0464, 0.1998, 0.0222, -1.0919, "
            "2.6851, 0.3740, 1.0000, 0.4739, 0.0114, 0.0243, 0.0069, 310.7211, 916.1947"
        )

        start_date = dt.datetime(2006, 1, 1)
        end_date = dt.datetime(2007, 12, 31)

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
            "elevation={elevation};"
            "rain_snow_fraction={rain_snow_fraction};"
            "nc_spec={tasmax};"
            "nc_spec={tasmin};"
            "nc_spec={pr}".format(
                ts=ts,
                params=params,
                start_date=start_date,
                end_date=end_date,
                init="155,455",
                name="Salmon",
                run_name="test-hmets-NRCAN",
                area="4250.6",
                elevation="843.0",
                latitude=lat,
                longitude=lon,
                rain_snow_fraction="RAINSNOW_DINGMAN",
                tasmax=json.dumps({"tasmax": {"linear_transform": (1.0, -273.15)}}),
                tasmin=json.dumps({"tasmin": {"linear_transform": (1.0, -273.15)}}),
                pr=json.dumps({"pr": {"linear_transform": (86400.0, 0.0)}}),
            )
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="raven-hmets",
            datainputs=datainputs,
        )

        assert_response_success(resp)

        # For testing purposes
        # out = get_output(resp.xml)
        # tmp_file, _ = urlretrieve(out['hydrograph'])
        # tmp=xr.open_dataset(tmp_file)
        # plt.plot(tmp['q_sim'])
        # plt.show()
