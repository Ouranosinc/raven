from pathlib import Path

from matplotlib import pyplot as plt
from pywps import FORMATS, ComplexInput, ComplexOutput, Format, LiteralInput, Process
from ravenpy.utilities.graphs import ts_graphs


class GraphIndicatorAnalysis(Process):
    def __init__(self):
        inputs = [
            ComplexInput(
                "ts_stats",
                "Stream flow time series statistics calculated by wps_tsstats",
                abstract="Stream flow time-series statistics",
                supported_formats=[FORMATS.NETCDF, Format(mime_type="application/zip")],
            ),
            LiteralInput(
                "trend",
                "Show trend in data using Thiel-Sen method",
                abstract="Display trend in time-series statistics",
                data_type="boolean",
                default=True,
                min_occurs=0,
                max_occurs=1,
            ),
            LiteralInput(
                "alpha",
                "test for trend using Mann-Kendall test",
                abstract="Test for trend in time-series statistics",
                data_type="float",
                default=0.05,
                min_occurs=0,
                max_occurs=1,
            ),
        ]

        outputs = [
            ComplexOutput(
                "graph_ts_stats",
                "Figure showing the time-series of the indices.",
                abstract="Time-series of statistical indices ",
                as_reference=True,
                supported_formats=(Format(mime_type="image/png"),),
            ),
        ]

        super().__init__(
            self._handler,
            identifier="ts_stats_graph",
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
        sim_fn = request.inputs["ts_stats"][0].file
        tr = request.inputs["trend"][0].data
        alph = request.inputs["alpha"][0].data

        # Create and save graphic
        fig = ts_graphs(sim_fn, trend=tr, alpha=alph)

        fig_ts_stats = Path(self.workdir) / "ts_graphs.png"
        fig.savefig(fig_ts_stats)
        plt.close(fig)

        response.outputs["graph_ts_stats"].file = str(fig_ts_stats)

        return response
