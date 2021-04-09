#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 13:23:00 2019

@author: ets
"""

from pathlib import Path

from matplotlib import pyplot as plt
from pywps import FORMATS, ComplexInput, ComplexOutput, Format, Process
from ravenpy.utilities.graphs import hydrograph, mean_annual_hydrograph


class GraphObjectiveFunctionFitProcess(Process):
    def __init__(self):
        inputs = [
            ComplexInput(
                "sims",
                "NetCDF containing q_sim and q_obs for model calibration fit check.",
                abstract="Stream flow simulation time series",
                supported_formats=[FORMATS.NETCDF],
            ),
        ]

        outputs = [
            ComplexOutput(
                "graph_objfun_fit",
                "Figure showing the observed and simulated streamflows",
                abstract="",
                as_reference=True,
                supported_formats=(Format(mime_type="image/png"),),
            ),
            ComplexOutput(
                "graph_objfun_annual_fit",
                "Figure showing the fit on the mean annual hydrograph.",
                abstract="",
                as_reference=True,
                supported_formats=(Format(mime_type="image/png"),),
            ),
        ]

        super(GraphObjectiveFunctionFitProcess, self).__init__(
            self._handler,
            identifier="graph_objective_function_fit",
            title="",
            version="1.0",
            abstract="",
            metadata=[],
            inputs=inputs,
            outputs=outputs,
            keywords=[],
            status_supported=True,
            store_supported=True,
        )

    def _handler(self, request, response):
        sim_fn = request.inputs["sims"][0].file

        # Create and save graphic
        fig = mean_annual_hydrograph([sim_fn])
        fig_fn_annual = Path(self.workdir) / "graph_objfun_annual_fit.png"
        fig.savefig(fig_fn_annual)
        plt.close(fig)

        fig = hydrograph([sim_fn])
        fig_fn_simple = Path(self.workdir) / "graph_objfun_fit.png"
        fig.savefig(fig_fn_simple)
        plt.close(fig)

        response.outputs["graph_objfun_fit"].file = str(fig_fn_simple)
        response.outputs["graph_objfun_annual_fit"].file = str(fig_fn_annual)

        return response
