import os
from pywps import Process
from pywps import LiteralInput
from pywps import ComplexInput, ComplexOutput
from pywps import FORMATS
# from pywps.app.Common import Metadata
import subprocess
from . import ravenio

import logging
from collections import OrderedDict as Odict
LOGGER = logging.getLogger("PYWPS")


class RavenProcess(Process):
    identifier = 'raven'
    abstract = 'Raven hydrological framework'
    title = ''
    version = ''

    nc = ComplexInput('nc', 'netCDF input files',
                      abstract='NetCDF file or files storing'
                               ' daily liquid precipitation (pr), '
                               'solid precipitation (prsn), '
                               'minimum temperature (tasmin), '
                               'maximum temperature (tasmax), '
                               'potential evapotranspiration (evspsbl) and '
                               'observed streamflow (qobs [m3/s]).',
                      min_occurs=1,
                      supported_formats=[FORMATS.NETCDF, FORMATS.TEXT])

    conf = ComplexInput('conf', 'Raven configuration files',
                        abstract="Raven rv files",
                        min_occurs=5,
                        max_occurs=5,
                        supported_formats=[FORMATS.TEXT])

    hydrograph = ComplexOutput('hydrograph', 'Hydrograph time series (mm)',
                               supported_formats=[FORMATS.NETCDF],
                               abstract='A netCDF file containing the outflow hydrographs (in m3/s) for all subbasins'
                                        'specified as `gauged` in the .rvh file. It reports period-ending time-'
                                        'averaged flows for the preceding time step, as is consistent with most '
                                        'measured stream gauge data (again, the initial flow conditions at the '
                                        'start of the first time step are included). If observed hydrographs are '
                                        'specified, they will be output adjacent to the corresponding modelled  '
                                        'hydrograph. ',
                               as_reference=True)

    storage = ComplexOutput('storage', 'Watershed storage time series (mm)',
                            abstract='A netCDF file describing the total storage of water (in mm) in all water '
                                     'storage compartments for each time step of the simulation. Mass balance '
                                     'errors, cumulative input (precipitation), and output (channel losses) are '
                                     'also included. Note that the precipitation rates in this file are '
                                     'period-ending, i.e., this is the precipitation rate for the time step '
                                     'preceding the time stamp; all water storage variables represent '
                                     'instantaneous reports of the storage at the time stamp indicate.',
                            supported_formats=[FORMATS.NETCDF],
                            as_reference=True)

    solution = ComplexOutput('solution', 'solution.rvc file to restart another simulation with the conditions '
                                         'at the end of this simulation.',
                             supported_formats=[FORMATS.TEXT],
                             as_reference=True)

    diagnostics = ComplexOutput('diagnostics', 'Performance diagnostic values',
                                abstract="Model diagnostic CSV file.",
                                supported_formats=[FORMATS.TEXT],
                                as_reference=True)
    def __init__(self):

        inputs = [self.nc, self.conf]

        outputs = [self.hydrograph, self.storage, self.solution, self.diagnostics]

        super(RavenProcess, self).__init__(
            self._handler,
            identifier=self.identifier,
            title=self.title,
            version=self.version,
            abstract=self.abstract,
            inputs=inputs,
            outputs=outputs,
            status_supported=True,
            store_supported=True
        )

    def _setup_dir(self):
        """Create directory structure to store model input files, executable and output results.

        workdir/  # Created by PyWPS. Is considered the model path.
           output/

        """
        output_path = os.path.join(self.workdir, 'output')

        # Create subdirectory
        os.mkdir(output_path)

    def _handler(self, request, response):

        self._setup_dir()
        config = {}
        files = [inpt.file for inpt in request.inputs['conf']]
        for f in files:
            name, ext = os.path.splitext(os.path.split(f)[-1])
            config[ext[1:]] = f

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

        out_files = {'hydrograph': 'Hydrographs.nc',
                     'storage': 'WatershedStorage.nc',
                     'solution': 'solution.rvc',
                     'diagnostics': 'Diagnostics.csv'}

        # Assign the response outputs to the full names
        for name in expected:
            fn = out_files[name]
            for f in files:
                if fn in f:
                    out[name] = os.path.join(out_dir, f)

            if name not in out:
                raise ValueError("No file named {} was found in output directory.".format(fn))

        return out


"""
    rvi = ComplexInput('rvi', 'Primary input file',
                       abstract="The primary input file stores the model simulation options and numerical options.",
                       min_occurs=1,
                       max_occurs=1,
                       supported_formats=[FORMATS.TEXT])

    rvp = ComplexInput('rvp', 'Classed parameter input file',
                       abstract="The classed parameter input file stores a database of soil, vegetation, river, "
                                "aquifer, and land class pro-perties. Not all classes specified in the *.rvp file "
                                "need to be included in the model.",
                       min_occurs=1,
                       max_occurs=1,
                       supported_formats=[FORMATS.TEXT])

    rvh = ComplexInput('rvh', 'HRU / Basin definition file',
                       abstract="The HRU/basin definition file describes the topology of the basin network and the "
                                "class membership of all constituent HRUs.",
                       min_occurs=1,
                       max_occurs=1,
                       supported_formats=[FORMATS.TEXT])

    rvt = ComplexInput('rvt', 'Time series input file',
                       abstract="The time series input file is used to store time series of forcing functions ("
                                "precipitation, temperature, etc.).",
                       min_occurs=1,
                       max_occurs=1,
                       supported_formats=[FORMATS.TEXT])

    rvc = ComplexInput('rvc', 'Initial conditions input file',
                       abstract="The initial conditions input file is used to store the initial conditions for the "
                                "model. By default, the initial conditions for all model state variables is zero, "
                                "and there are no required commands in this file (it could even be completely "
                                "empty).",
                       min_occurs=0,
                       default="",
                       supported_formats=[FORMATS.TEXT])
"""
