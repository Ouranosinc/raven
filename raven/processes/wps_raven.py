import logging
from collections import defaultdict
from pathlib import Path
import json

from pywps import Process, Format, LiteralOutput
from raven.utils import archive_sniffer, single_file_check
from raven.models import Raven
from . import wpsio as wio

LOGGER = logging.getLogger("PYWPS")


class RavenProcess(Process):
    identifier = 'raven'
    abstract = 'Raven hydrological framework'
    title = "Run the Raven hydrological framework using model configuration files and forcing time series. In " \
            "the `rvt` file, only provide the name of the forcing file, not an absolute or relative path."
    version = '0.1'

    tuple_inputs = {}
    inputs = [wio.ts, wio.nc_spec, wio.conf]
    outputs = [wio.hydrograph, wio.storage, wio.solution, wio.diagnostics, wio.rv_config]
    model_cls = Raven

    def __init__(self):

        super(RavenProcess, self).__init__(
            self._handler,
            identifier=self.identifier,
            title=self.title,
            version=self.version,
            abstract=self.abstract,
            inputs=self.inputs,
            outputs=self.outputs,
            status_supported=True,
            store_supported=True
        )

    def model(self, request):
        """Return model class."""
        return self.model_cls(workdir=self.workdir)

    def meteo(self, request):
        """Return meteo input files."""
        return [f.file for f in request.inputs.pop('ts')]

    def region(self, request):
        """Return region shape file."""
        extensions = ['.gml', '.shp', '.gpkg', '.geojson', '.json']
        region_vector = request.inputs.pop('region_vector')[0].file
        return single_file_check(archive_sniffer(region_vector, working_dir=self.workdir, extensions=extensions))

    def options(self, request):
        """Parse model options."""
        # Input specs dictionary. Could be all given in the same dict or a list of dicts.
        kwds = defaultdict(list)
        for spec in request.inputs.pop('nc_spec', []):
            kwds.update(json.loads(spec.data))

        # Parse all other input parameters
        for name, objs in request.inputs.items():
            for obj in objs:

                # Namedtuples
                if name in self.tuple_inputs:
                    data =self.parse_tuple(obj)

                # Other parameters
                else:
                    data = obj.data

                if name in Raven._parallel_parameters:
                    kwds[name].append(data)
                else:
                    kwds[name] = data

        return kwds

    def parse_tuple(self, obj):
        csv = obj.data.replace('(', '').replace(')', '')
        arr = map(float, csv.split(','))
        return self.tuple_inputs[obj.identifier](*arr)

    def run(self, model, ts, kwds):
        """Run the model.

        If keywords contain `rvc`, initialize the model using the initial condition file."""
        model(ts=ts, **kwds)

    def initialize(self, model, request):
        """Set initial conditions from a solution.rvc file.

        This is used by emulators.
        """
        rvc = self.get_config(request, ids=("rvc",))["rvc"]
        if rvc:
            model.resume(rvc)

    def _handler(self, request, response):
        response.update_status('PyWPS process {} started.'.format(self.identifier), 0)

        model = self.model(request)

        # Model configuration (RV files)
        config = self.get_config(request, ids=("conf",))
        if config:
            if len(config) > 1:
                raise NotImplementedError("Multi-model simulations are not yet supported.")
            conf = list(config.values()).pop()
            model.configure(conf.values())

        self.initialize(model, request)

        # Input data files
        ts = self.meteo(request)

        # Model options
        kwds = self.options(request)

        # Launch model with input files
        self.run(model, ts, kwds)

        # Store output files name. If an output counts multiple files, they'll be zipped.
        for key in response.outputs.keys():
            val = model.outputs.get(key)
            if val is not None:
                if isinstance(response.outputs[key], LiteralOutput):
                    response.outputs[key].data = str(val)
                else:
                    response.outputs[key].file = str(val)
                    if val.suffix == '.zip':
                        response.outputs[key].data_format = \
                            Format('application/zip', extension='.zip', encoding='base64')

        return response

    def get_config(self, request, ids=("conf",)):
        """Return a dictionary storing the configuration files content."""
        config = defaultdict(dict)
        for key in ids:
            if key in request.inputs:
                conf = request.inputs.pop(key)
                for obj in conf:
                    fn = Path(obj.file)
                    config[fn.stem][fn.suffix[1:]] = fn

        return config
