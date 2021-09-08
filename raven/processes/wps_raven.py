import json
import logging
import string
import traceback
from collections import defaultdict
from pathlib import Path

from pywps import Format, LiteralOutput, Process
from pywps.app.exceptions import ProcessError
from ravengis.io import archive_sniffer
from ravenpy.config.commands import HRU
from ravenpy.models import Raven

from raven.utilities import single_file_check

from . import wpsio as wio

LOGGER = logging.getLogger("PYWPS")


class RavenProcess(Process):
    identifier = "raven"
    abstract = "Raven hydrological framework"
    title = (
        "Run the Raven hydrological framework using model configuration files and forcing time series. In "
        "the `rvt` file, only provide the name of the forcing file, not an absolute or relative path."
    )
    version = "0.1"

    tuple_inputs = dict()
    inputs = [wio.ts, wio.nc_spec, wio.conf]
    outputs = [
        wio.hydrograph,
        wio.storage,
        wio.solution,
        wio.diagnostics,
        wio.rv_config,
    ]
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
            store_supported=True,
        )

    def model(self, request):
        """Return model class."""
        return self.model_cls(workdir=self.workdir)

    def meteo(self, request):
        """Return meteo input files."""
        return [f.file for f in request.inputs.pop("ts")]

    def region(self, request):
        """Return region shape file."""
        extensions = [".gml", ".shp", ".gpkg", ".geojson", ".json"]
        region_vector = request.inputs.pop("region_vector")[0].file
        return single_file_check(
            archive_sniffer(
                region_vector, working_dir=self.workdir, extensions=extensions
            )
        )

    def options(self, request):
        """Parse model options."""
        # Input specs dictionary. Could be all given in the same dict or a list of dicts.
        kwds = defaultdict(list)
        for spec in request.inputs.pop("nc_spec", []):
            kwds.update(json.loads(spec.data))

        # This is an alternate interface to the legacy one, composed of: area/latitude/longitude/elevation
        for hrus in request.inputs.pop("hrus", []):
            kwds["hrus"] = [HRU(**attrs) for attrs in json.loads(hrus.data)]

        # Parse all other input parameters
        for name, objs in request.inputs.items():
            for obj in objs:

                # Namedtuple objects
                if name in self.tuple_inputs:
                    data = self.parse_tuple(obj)

                # Other parameters
                else:
                    data = obj.data

                if name in Raven._parallel_parameters:
                    kwds[name].append(data)
                else:
                    kwds[name] = data

        return kwds

    def parse_tuple(self, obj):
        csv = obj.data.replace("(", "").replace(")", "")
        arr = map(float, csv.split(","))
        return self.tuple_inputs[obj.identifier](*arr)

    def run(self, model, ts, kwds):
        """Run the model.

        If keywords contain `rvc`, initialize the model using the initial condition file."""
        model(ts=ts, **kwds)

    def _handler(self, request, response):
        response.update_status(f"PyWPS process {self.identifier} started.", 0)

        model = self.model(request)

        # Model configuration (zipped RV files in `conf` input)
        if "conf" in request.inputs:
            model.configure(self.get_config(request).values())

        # Initial conditions (`rvc` input)
        if "rvc" in request.inputs:
            model.resume(request.inputs.pop("rvc")[0].file)

        if "random_numbers" in request.inputs:
            model.config.set_rv_file(request.inputs.pop("random_numbers")[0].file)

        # Input data files
        ts = self.meteo(request)

        # Model options
        kwds = self.options(request)

        # Launch model with input files
        try:
            self.run(model, ts, kwds)
        except Exception as exc:
            LOGGER.exception(exc)
            err_msg = traceback.format_exc()
            # By default the error message is limited to 300 chars and strips
            # many special characters
            raise ProcessError(
                err_msg, max_length=len(err_msg), allowed_chars=string.printable
            ) from exc

        # Store output files name. If an output counts multiple files, they'll be zipped.
        for key in response.outputs.keys():
            val = model.outputs.get(key)
            if val is not None:
                if isinstance(response.outputs[key], LiteralOutput):
                    response.outputs[key].data = str(val)
                else:
                    response.outputs[key].file = str(val)
                    if val.suffix == ".zip":
                        response.outputs[key].data_format = Format(
                            "application/zip", extension=".zip", encoding="base64"
                        )
            else:
                response.outputs[key].data = ""

        return response

    def get_config(self, request):
        """Return a dictionary storing the configuration files content."""
        import zipfile

        config = defaultdict(dict)
        conf = request.inputs.pop("conf")[0].file

        z = zipfile.ZipFile(conf)
        z.extractall(self.workdir)

        for f in z.namelist():
            fn = Path(self.workdir) / f
            config[fn.stem][fn.suffix[1:]] = fn

        if len(config) > 1:
            raise NotImplementedError("Multi-model simulations are not yet supported.")
        return list(config.values()).pop()
