import zipfile
from pathlib import Path

from matplotlib import pyplot as plt
from pywps import ComplexInput, ComplexOutput
from pywps import FORMATS
from pywps import Format
from pywps import Process

from raven.utilities.graphs import mean_annual_hydrograph, hydrograph


class GraphEnsUncertaintyProcess(Process):
    def __init__(self):
        inputs = [ComplexInput('sims', 'Stream flow simulations ensemble',
                               abstract='Stream flow simulation time series',
                               supported_formats=[FORMATS.NETCDF, Format(mime_type='application/zip')]),
                  ]

        outputs = [
            ComplexOutput('graph_ensemble_hydrographs', 'Figure showing the simple hydrographs of the included models.',
                          abstract="",
                          as_reference=True,
                          supported_formats=(Format(mime_type='image/png'),)),

            ComplexOutput('graph_annual_hydrographs', 'Figure showing the spread for the mean annual hydrograph.',
                          abstract="",
                          as_reference=True,
                          supported_formats=(Format(mime_type='image/png'),)),
        ]

        super(GraphEnsUncertaintyProcess, self).__init__(
            self._handler,
            identifier="graph_ensemble_uncertainty",
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
        sim_fn = request.inputs['sims'][0].file

        # Extract files from archive in temp directory
        tmp = Path(self.workdir) / 'sims'
        with zipfile.ZipFile(sim_fn) as z:
            z.extractall(tmp)

        # Create and save graphic
        fig = mean_annual_hydrograph(tmp.glob('*.nc'))
        fig_fn_annual = Path(self.workdir) / 'ensemble_annual_hydrographs.png'
        fig.savefig(fig_fn_annual)
        plt.close(fig)

        fig = hydrograph(tmp.glob('*.nc'))
        fig_fn_simple = Path(self.workdir) / 'simple_hydrographs.png'
        fig.savefig(fig_fn_simple)
        plt.close(fig)

        response.outputs['graph_ensemble_hydrographs'].file = str(fig_fn_simple)
        response.outputs['graph_annual_hydrographs'].file = str(fig_fn_annual)

        return response
