from pywps import Process
from raven.models import Raven
from . import wpsio as wio
import logging
from pathlib import Path
from collections import OrderedDict as Odict, defaultdict
LOGGER = logging.getLogger("PYWPS")


class RavenProcess(Process):
    identifier = 'raven'
    abstract = 'Raven hydrological framework'
    title = "Run the Raven hydrological framework using model configuration files and forcing time series. In " \
            "the `rvt` file, only provide the name of the forcing file, not an absolute or relative path."
    version = '0.1'

    tuple_inputs = {}
    inputs = [wio.ts, wio.conf]
    outputs = [wio.hydrograph, wio.storage, wio.solution, wio.diagnostics]
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

    def _handler(self, request, response):
        response.update_status('PyWPS process {} started.'.format(self.identifier), 0)

        model = self.model_cls(self.workdir)

        # Model configuration
        if 'conf' in request.inputs:
            conf = request.inputs.pop('conf')
            config = self.get_config(conf)
            model.configure(config.values())

        # Input data files
        ts = [f.file for f in request.inputs.pop('ts')]

        # Parse all other input parameters
        for name, obj in request.inputs.items():

            # Namedtuples
            if name in self.tuple_inputs:
                arr = map(float, obj[0].data.split(','))
                data = self.tuple_inputs[name](*arr)

            # Other parameters
            else:
                data = obj[0].data

            model.assign(name, data)

        # Launch model with input files
        model.run(ts=ts)

        for key, val in model.outputs.items():
            response.outputs[key].file = str(val)

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
