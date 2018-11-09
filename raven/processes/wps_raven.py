import os
from pywps import Process
import subprocess
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

    def _setup_dir(self, name):
        """Create directory structure to store model input files, executable and output results.

        Model configuration files and time series inputs are stored directly in the working directory.

        workdir/  # Created by PyWPS. Is considered the model path.
           output/

        """
        output_path = os.path.join(self.workdir, name, 'output')
        model_path = os.path.join(self.workdir, name, 'model')

        # Create subdirectory
        os.mkdir(output_path)
        os.mkdir(model_path)

    def _handler(self, request, response):
        response.update_status('PyWPS process {self.identifier} started.'.format(self), 0)

        # Create output directory.
        self._setup_dir()

        # Read configuration files
        config = defaultdict(dict)
        for input in request.inputs['conf']:
            file = Path(input.file)
            config[file.stem][file.suffix[1:]] = file




        cmd = os.path.join(self.workdir, 'raven')
        os.symlink(ravenio.raven_exec, cmd)

        # Run the simulation
        subprocess.call([cmd, os.path.join(self.workdir, name), '-o', os.path.join(self.workdir, 'output')])
        cmdline = " ".join([cmd, os.path.join(self.workdir, name), '-o', os.path.join(self.workdir, 'output')])
        # Get the output files.
        outs = self._match_outputs(response.outputs.keys())
        for key, val in outs.items():
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

