import logging
from collections import defaultdict
from pathlib import Path
import json

from pywps import Process, Format, LiteralOutput

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
        return self.model_cls(workdir=self.workdir)

    def _handler(self, request, response):
        response.update_status('PyWPS process {} started.'.format(self.identifier), 0)

        model = self.model(request)

        # Model configuration
        if 'conf' in request.inputs:
            conf = request.inputs.pop('conf')
            config = self.get_config(conf)
            model.configure(config.values())

        # Input data files
        ts = [f.file for f in request.inputs.pop('ts')]

        # Input specs dictionary. Could a all given in the same dict or a list of dicts.
        nc_spec = {}
        for spec in request.inputs.pop('nc_spec', []):
            nc_spec.update(json.loads(spec.data))

        for key, val in nc_spec.items():
            model.assign(key, val)

        # Parse all other input parameters
        kwds = defaultdict(list)
        for name, objs in request.inputs.items():
            for obj in objs:

                # Namedtuples
                if name in self.tuple_inputs:
                    csv = obj.data.replace('(', '').replace(')', '')
                    arr = map(float, csv.split(','))
                    data = self.tuple_inputs[name](*arr)

                # Other parameters
                else:
                    data = obj.data

                if name in model._parallel_parameters:
                    kwds[name].append(data)
                else:
                    model.assign(name, data)

        # Launch model with input files
        model(ts=ts, **kwds)

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

    @staticmethod
    def get_config(conf):
        """Return a dictionary storing the configuration files content."""
        config = defaultdict(dict)
        for obj in conf:
            fn = Path(obj.file)
            config[fn.stem][fn.suffix[1:]] = fn

        if len(config.keys()) > 1:
            raise NotImplementedError("Multi-model simulations are not yet supported.")

        return config[fn.stem]
