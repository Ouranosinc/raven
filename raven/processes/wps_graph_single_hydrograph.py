#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 14:31:29 2019

@author: ets
"""

from pathlib import Path

from matplotlib import pyplot as plt
from pywps import ComplexInput, ComplexOutput
from pywps import FORMATS
from pywps import Format
from pywps import Process

from raven.utilities.graphs import mean_annual_hydrograph, hydrograph, spaghetti_annual_hydrograph


class GraphSingleHydrographProcess(Process):
    def __init__(self):
        inputs = [ComplexInput('sim', 'Stream flow simulations NetCDF File',
                               abstract='Stream flow simulation time series',
                               supported_formats=[FORMATS.NETCDF]),
                  ]

        outputs = [
            ComplexOutput('graph_single_hydrographs', 'Figure showing the simple hydrographs of the included models.',
                          abstract="",
                          as_reference=True,
                          supported_formats=(Format(mime_type='image/png'),)),

            ComplexOutput('graph_annual_hydrographs', 'Figure showing the spread for the mean annual hydrograph.',
                          abstract="",
                          as_reference=True,
                          supported_formats=(Format(mime_type='image/png'),)),

            ComplexOutput('graph_spaghetti_hydrographs', 'Figure showing the annual hydrographs for each year.',
                          abstract="",
                          as_reference=True,
                          supported_formats=(Format(mime_type='image/png'),)),
        ]

        super(GraphSingleHydrographProcess, self).__init__(
            self._handler,
            identifier="graph_single_hydrograph",
            title="",
            version="1.0",
            abstract="",
            metadata=[],
            inputs=inputs,
            outputs=outputs,
            keywords=[],
            status_supported=True,
            store_supported=True)

    def _handler(self, request, response):
        sim_fn = Path(request.inputs['sim'][0].workdir).glob('*.nc')

        # Create and save graphics
        fig = mean_annual_hydrograph(sim_fn)
        fig_fn_annual = Path(self.workdir) / 'single_annual_hydrographs.png'
        fig.savefig(fig_fn_annual)
        plt.close(fig)

        sim_fn = Path(request.inputs['sim'][0].workdir).glob('*.nc')

        fig = hydrograph(sim_fn)
        fig_fn_simple = Path(self.workdir) / 'simple_hydrograph.png'
        fig.savefig(fig_fn_simple)
        plt.close(fig)

        sim_fn = request.inputs['sim'][0].file

        fig = spaghetti_annual_hydrograph(sim_fn)
        fig_fn_spag = Path(self.workdir) / 'spaghetti_hydrographs.png'
        fig.savefig(fig_fn_spag)
        plt.close(fig)

        response.outputs['graph_single_hydrographs'].file = str(fig_fn_simple)
        response.outputs['graph_annual_hydrographs'].file = str(fig_fn_annual)
        response.outputs['graph_spaghetti_hydrographs'].file = str(fig_fn_spag)

        return response
