import datetime as dt
import os
import tempfile
import numpy as np
import pandas as pd
import xarray as xr
import pytest
import matplotlib.pyplot as plt
from copy import deepcopy
from raven.utilities.data_assimilation import assimilate, perturbation
from raven.models import GR4JCN
from raven.models.rv import RVC
from raven.models.state import BasinStateVariables
from .common import TESTDATA
import pdb


def test_perturbation():
    ts = TESTDATA["raven-gr4j-cemaneige-nc-ts"]
    ds = xr.open_dataset(ts)

    tmax = ds.tmax.isel(time=slice(0, 10))
    p_tmax = perturbation(tmax, "norm", 0.01, members=50)
    np.testing.assert_allclose(p_tmax.mean("members"), tmax, rtol=0.1)

    rain = ds.rain.isel(time=slice(30, 60))
    p_rain = perturbation(rain, "gamma", 0.01, members=50)
    np.testing.assert_allclose(p_rain.mean("members"), rain, rtol=0.1)

    assert p_tmax.attrs == ds.tmax.attrs
    assert p_rain.attrs == ds.rain.attrs


class TestAssimilationGR4JCN:
    def test_simple(self):

        # get timeseries
        ts = TESTDATA["raven-gr4j-cemaneige-nc-ts"]

        # set number of members. Using 7 here to make it easier to find and debug.
        n_members = 7

        # Perturbation parameters for the assimilation, keyed by standard_name
        std = {
            "rainfall": 0.30,
            "prsn": 0.30,
            "tasmin": 2.0,
            "tasmax": 2.0,
            "water_volume_transport_in_river_channel": 0.15,
        }

        # Use the same random seed for both tasmin and tasmax
        rs = np.random.SeedSequence(None).generate_state(1)[0]
        seed = {"tasmin": rs,
                "tasmax": rs}

        # Perturbation distribution
        dists = {
            "pr": "gamma",
            "rainfall": "gamma",
            "prsn": "gamma",
            "water_volume_transport_in_river_channel": "rnorm"}

        qkey = "water_volume_transport_in_river_channel"
        if qkey not in std:
            raise ValueError("Assimilation requires perturbing the flow variable.")

        # Assimilation variables (from HRUStateVariable)
        assim_var = ("soil0", "soil1")

        # Assimilation periods
        assim_days = [10] + 4 * [3]

        # GR4JCN model instance
        model = GR4JCN()

        # set the start and end dates for the first assimilation period, warm-up
        start_date = dt.datetime(2000, 6, 1)
        end_date = start_date + dt.timedelta(days=sum(assim_days))

        # Set model options
        model.rvh.name = "Salmon"
        model.rvh.area = "4250.6"
        model.rvh.elevation = "843.0"
        model.rvh.latitude = 54.4848
        model.rvh.longitude = -123.3659

        model.rvp.params = model.params(
            0.1353389, -0.005067198, 576.8007, 6.986121, 1.102917, 0.9224778
        )  # SALMON

        # ==== Initialization (just to get reasonable states) ====
        # Set initialization run options
        model.rvi.run_name = "init"
        model.rvi.start_date = start_date
        # Richard: what is the end date policy for the init run ?
        model.rvi.end_date = start_date + dt.timedelta(days=assim_days[0])

        # Run the model
        model([ts])

        # Extract final model states
        hru_state, basin_state = model.get_final_state()
        xa = n_members * [getattr(hru_state, key) for key in assim_var]
        hru_states = n_members * [hru_state]
        basin_states = n_members * [basin_state]

        # === Create perturbed time series for full assimilation period ====
        perturbed = {}
        for key, s in std.items():
            nc = model.rvt.get(key)

            with xr.open_dataset(nc.path) as ds:
                da = ds.get(nc.var_name).sel(time=slice(start_date, end_date))

                perturbed[key] = perturbation(
                    da, dists.get(key, "norm"), std=s, seed=seed.get(key, None), member=n_members
                )

                # Save flow for later
                if key == qkey:
                    q_obs = da

        # Write to disk
        p_fn = model.workdir / "perturbed_forcing.nc"
        perturbed = xr.Dataset(perturbed)
        perturbed.to_netcdf(p_fn, mode="w")

        # ==== Assimilation ====
        q_assim = []
        sd = start_date
        for i, ndays in enumerate(assim_days):

            dates = [sd + dt.timedelta(days=x) for x in range(ndays)]
            model.rvi.end_date = dates[-1]
            model.rvi.run_name = f"assim_{i}"

            # Perform the first assimilation step here
            [xa, model] = assimilate(model, p_fn, q_obs, assim_var, basin_states, hru_states, dates)

            # Save streamflow simulation
            q_assim.append(model.q_sim.isel(nbasins=0))

            # Update the start-time for the next loop
            sd += dt.timedelta(days=ndays)
            model.rvi.start_date = sd

            # Get new initial conditions and feed assimilated values
            hru_states, basin_states = model.get_final_state()
            hru_states = [hru_states[i]._replace(**dict(zip(assim_var, xa[:, i]))) for i in range(n_members)]

        q_assim = xr.concat(q_assim, dim="time")

        # ==== Reference run ====
        model.rvi.run_name = "ref"
        model.rvi.start_date = start_date
        model.rvi.end_date = end_date
        model.rvc = RVC(soil0=None, soil1=15, basin_state=BasinStateVariables())
        model([ts, ])

        # We can now plot everything!
        plt.plot(q_assim.T, "r", label="Assimilated")  # plot the assimilated flows
        plt.plot(q_obs.T, "b", label="Observed")  # plot the observed flows
        plt.plot(
            model.q_sim, "g", label="Simulated"
        )  # plot the open_loop (simulation with no assimilation)
        # plt.legend()
        plt.show()

        assert q_assim.shape[0] == n_members
        assert q_assim.shape[1] == sum(assim_days)
