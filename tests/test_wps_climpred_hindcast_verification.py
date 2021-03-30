import datetime as dt
from urllib.request import urlretrieve

from pywps import Service
from pywps.tests import assert_response_success
from ravenpy.utilities.testdata import get_local_testdata
from ravenpy.utilities.forecasting import make_ESP_hindcast_dataset
from raven.processes import ClimpredHindcastVerificationProcess
import tempfile
from .common import CFG_FILE, client_for, get_output


class TestHindcastEvaluationProcess:
    def test_hindcast_evaluation_rank_histo(self):
        client = client_for(
            Service(
                processes=[ClimpredHindcastVerificationProcess()], cfgfiles=CFG_FILE
            )
        )

        # Prepare the model parameters and forecast details
        model = "GR4JCN"
        params = (0.529, -3.396, 407.29, 1.072, 16.9, 0.947)
        forecast_duration = 3
        ts = get_local_testdata(
            "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
        )
        rvc = get_local_testdata("gr4j_cemaneige/solution.rvc")
        # Make the hindcasts for each initialization date. Here we will extract
        # ESP forecasts for a given calendar date for the years in "included_years"
        # as hindcast dates. Each ESP hindcast uses all available data in the ts dataset,
        # so in this case we will have 56/57 members for each hindcast initialization
        # depending on the date that we start on. The "hindcasts" dataset contains
        # all of the flow data from the ESP hindcasts for the initialization dates.
        # The "qobs" dataset contains all qobs in the timeseries: Climpred will
        # sort it all out during its processing. Note that the format of these datasets
        # is tailor-made to be used in climpred, and thus has specific dimension names.
        hcasts, obs = make_ESP_hindcast_dataset(
            model_name=model,
            forecast_date=dt.datetime(1955, 6, 30),
            included_years=list(range(1957, 1959)),
            forecast_duration=forecast_duration,
            ts=ts,
            area="4250.6",
            elevation="843.0",
            latitude=54.4848,
            longitude=-123.3659,
            params=params,
            rvc=str(rvc),
        )

        tmpdir = tempfile.mkdtemp()
        tmpfile_hcst = tmpdir + "/hcst_wps.nc"
        tmpfile_obs = tmpdir + "/qobs_wps.nc"

        hcasts.to_netcdf(tmpfile_hcst)
        obs.to_netcdf(tmpfile_obs)
        metric = "rank_histogram"

        datainputs = (
            "hindcasts=files@xlink:href=file://{hindcasts};"
            "observations=files@xlink:href=file://{observations};"
            "metric={metric}".format(
                hindcasts=tmpfile_hcst,
                observations=tmpfile_obs,
                metric=metric,
            )
        )

        resp = client.get(
            service="WPS",
            request="Execute",
            version="1.0.0",
            identifier="climpred_hindcast_verification",
            datainputs=datainputs,
        )

        assert_response_success(resp)
        out = get_output(resp.xml)
        metrics, _ = urlretrieve(out["verification_metrics"])
