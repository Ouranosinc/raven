from pywps import ComplexInput, ComplexOutput
from pywps import Format
from pywps import Process
import zipfile
from pathlib import Path
from raven.utilities.graphs import ensemble_uncertainty_annual
from matplotlib import pyplot as plt
import os

class GraphEnsUncertaintyProcess(Process):
    def __init__(self):
        inputs = [ComplexInput('sims', 'Stream flow simulations ensemble',
                               abstract='Stream flow simulation time series',
                               supported_formats=(Format(mime_type='application/zip'),)),
                  ]

        outputs = [ComplexOutput('graphic', 'Figure showing the spread for the mean annual hydrograph.',
                                 abstract="",
                                 as_reference=True,
                                 supported_formats=(Format(mime_type='image/png'), )),
                   ]

        super(GraphEnsUncertaintyProcess, self).__init__(
            self._handler,
            identifier="graph-ensemble-uncertainty",
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
        nameList=os.listdir(tmp) # Get the filenames from the tmp folder
        fig = ensemble_uncertainty_annual(tmp.glob('*.nc'),nameList)
        fig_fn = Path(self.workdir) / 'ensemble_uncertainty.png'
        fig.savefig(fig_fn)
        plt.close(fig)

        response.outputs['graphic'].file = str(fig_fn)
        return response
