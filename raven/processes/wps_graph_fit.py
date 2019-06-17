
from pathlib import Path

from matplotlib import pyplot as plt
# from plotly.tools import mpl_to_plotly
from pywps import ComplexInput, LiteralInput, ComplexOutput
from pywps import FORMATS
from pywps import Format
from pywps import Process
import xarray as xr

from raven.utilities.graphs import ts_fit_graph


class GraphFitProcess(Process):
    def __init__(self):
        inputs = [ComplexInput('ts', 'Stream flow time series',
                               abstract='Stream flow time series',
                               supported_formats=[FORMATS.NETCDF]),
                  ComplexInput('params', 'Distribution parameters',
                               abstract='Statistical distribution parameters fitted to time series',
                               supported_formats=[FORMATS.NETCDF]),
                  LiteralInput('variable', "Variable name",
                               abstract="Name of time series variable. If none will default to the first data variable"
                                        "found in file.",
                               data_type='string',
                               min_occurs=0,
                               default=""),
                  LiteralInput('format', "Output graphic format",
                               abstract="Graphic format.",
                               data_type='string',
                               default='png',
                               min_occurs=0,
                               allowed_values=['png', 'jpeg', 'pdf'])
                  ]

        outputs = [
            ComplexOutput('graph_fit', 'Graphic',
                          abstract="Graphic showing time series histogram and the probability density "
                                   "function of the fitted distribution.",
                          as_reference=True,
                          supported_formats=(Format(mime_type='image/png'),
                                             Format(mime_type='image/jpeg'),
                                             Format(mime_type='application/pdf'),
                                             Format(mime_type='application/json'),
                                             )),
        ]

        super(GraphFitProcess, self).__init__(
            self._handler,
            identifier="graph_fit",
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
        ts_fn = request.inputs['ts'][0].file
        p_fn = request.inputs['params'][0].file
        v = request.inputs['variable'][0].data
        format = request.inputs['format'][0].data

        # Create and save graphics
        ds = xr.open_dataset(ts_fn)
        if v == '':
            v = list(ds.data_vars.keys())[0]
        ts = ds[v]

        p = xr.open_dataset(p_fn)['params']  # Name of variable is hard-coded

        fig = ts_fit_graph(ts, p)

        if format == 'plotly':
            # This is not working great with this figure due to the twin axes.
            raise NotImplementedError
            # Create plotly object
            # obj = mpl_to_plotly(fig)

            # Convert to JSON
            # response.outputs['graph_fit'].data = obj.to_json()
            # response.outputs['graph_fit'].data_format = Format('application/json')

        else:
            fig_fn = Path(self.workdir) / ('ts_fit.' + format)
            fig.savefig(fig_fn, format=format)
            plt.close(fig)
            response.outputs['graph_fit'].file = str(fig_fn)
            if format in ['png', 'jpeg']:
                response.outputs['graph_fit'].data_format = Format('image/{}'.format(format))
            elif format == 'pdf':
                response.outputs['graph_fit'].data_format = Format('application/pdf')

        return response
