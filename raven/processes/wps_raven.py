import os
from pywps import Process
import subprocess
from raven.models import Raven
from . import ravenio
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

    inputs = [wio.ts, wio.conf]
    outputs = [wio.hydrograph, wio.storage, wio.solution, wio.diagnostics]

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

        model = Raven(self.workdir)

        # Read configuration files
        config = defaultdict(dict)
        for obj in request.inputs['conf']:
            fn = Path(obj.file)
            config[fn.stem][fn.suffix[1:]] = fn

        if len(config.keys()) > 1:
            raise NotImplementedError("Multi-model simulations are not yet supported.")

        model.configure(**config[fn.stem])
        ts = [f.file for f in request.inputs['ts']]
        model.run(ts=ts)

        for key, val in model.outputs.items():
            response.outputs[key].file = val

        return response

    def _match_outputs(self, expected):
        """Match actual output files to known expected files.

        Return a dictionary of file paths for each expected input.
        """
        import glob

        out = {}
        out_dir = os.path.join(self.workdir, 'output')
        files = glob.glob(os.path.join(out_dir, '*'))

        # Assign the response outputs to the full names
        for name in expected:
            fn = ravenio.output_filenames[name]
            for f in files:
                if fn in f:
                    out[name] = os.path.join(out_dir, f)

            if name not in out:
                raise ValueError("No file named {} was found in output directory.".format(fn))

        return out

