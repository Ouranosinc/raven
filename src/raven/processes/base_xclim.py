import logging
import os
from functools import reduce
from itertools import cycle
from operator import mul

import requests
import xarray as xr
from anyascii import anyascii
from pywps import FORMATS, ComplexInput, ComplexOutput, LiteralInput, Process
from pywps.app.Common import Metadata

LOGGER = logging.getLogger("PYWPS")


def make_xclim_indicator_process(name, xci):
    """Create a WPS Process subclass from an xclim `Indicator` class instance."""
    attrs = xci.json()

    # Sanitize name
    # name = attrs['identifier'].replace('{', '_').replace('}', '_').replace('__', '_')

    process_class = type(
        str(name) + "Process",
        (_XclimIndicatorProcess,),
        {"xci": xci, "__doc__": attrs["abstract"]},
    )

    return process_class


class _XclimIndicatorProcess(Process):
    """Dummy xclim indicator process class.

    Set xci to the xclim indicator in order to have a working class"""

    xci = None

    def __init__(self):
        """Create a WPS process from an xclim indicator class instance."""

        if self.xci is None:
            raise AttributeError(
                "Use the `make_xclim_indicator_process` function instead."
            )

        attrs = self.xci.json()

        outputs = [
            ComplexOutput(
                "output",
                "Function output in netCDF",
                abstract="The indicator values computed on the original input grid.",
                as_reference=True,
                supported_formats=[FORMATS.NETCDF],
            ),
            ComplexOutput(
                "output_log",
                "Logging information",
                abstract="Collected logs during process run.",
                as_reference=True,
                supported_formats=[FORMATS.TEXT],
            ),
        ]

        identifier = attrs["identifier"]
        super().__init__(
            self._handler,
            identifier=identifier,
            version="0.1",
            title=anyascii(attrs["long_name"]),
            abstract=anyascii(attrs["abstract"]),
            inputs=self.load_inputs(eval(attrs["parameters"])),
            outputs=outputs,
            status_supported=True,
            store_supported=True,
        )

    def load_inputs(self, params):
        # Ideally this would be based on the Parameters docstring section rather than name conventions.
        inputs = []

        for name, attrs in params.items():
            if name in ["tas", "tasmin", "tasmax", "pr", "prsn"]:
                inputs.append(make_nc_input(name))
            if name in ["da", "q", "arr"]:
                inputs.append(make_nc_input(name))
                inputs.append(make_variable())
            elif name in ["tn10", "tn90", "t10", "t90"]:
                inputs.append(make_nc_input(name))
            elif name in ["thresh_tasmin", "thresh_tasmax"]:
                inputs.append(make_thresh(name, attrs["default"], attrs["desc"]))
            elif name in [
                "thresh",
            ]:
                inputs.append(make_thresh(name, attrs["default"], attrs["desc"]))
            elif name in [
                "freq",
            ]:
                inputs.append(make_freq(name, attrs["default"]))
            elif name in [
                "window",
            ]:
                inputs.append(make_window(name, attrs["default"], attrs["desc"]))
            elif name in [
                "mode",
            ]:
                inputs.append(make_mode(name, attrs["desc"]))
            elif name in [
                "op",
            ]:
                inputs.append(make_op(name, attrs["desc"]))
            elif name in [
                "t",
            ]:
                inputs.append(make_t(name, attrs["desc"]))
            elif name in [
                "dist",
            ]:
                inputs.append(make_dist(name, attrs["desc"]))
            elif name in [
                "indexer",
            ]:
                inputs.extend(make_indexer())
            else:
                # raise NotImplementedError(name)
                LOGGER.warning(f"not implemented: {name}")

        return inputs

    def try_opendap(self, request):
        """Try to open the file as an OPeNDAP url and chunk it"""
        url = request.url
        if url and not url.startswith("file"):
            r = requests.get(url + ".dds")
            if r.status_code == 200 and r.content.decode().startswith("Dataset"):
                ds = xr.open_dataset(url)
                chunks = chunk_dataset(ds, max_size=1000000)
                ds = ds.chunk(chunks)
                self.write_log(f"Opened dataset as an OPeNDAP url: {url}")
                return ds

        self.write_log(f"Downloading dataset for url: {url}")

        # accessing the file property loads the data in the data property
        # and writes it to disk
        filename = request.file
        # we need to clean up the data property
        # if we don't do this, it will be written in the database and
        # to the output status xml file and it can get too large
        request._data = ""

        return xr.open_dataset(filename)

    def log_file_path(self):
        return os.path.join(self.workdir, "log.txt")

    def write_log(self, message):
        open(self.log_file_path(), "a").write(message + "\n")
        LOGGER.info(message)

    def _handler(self, request, response):
        response.outputs["output_log"].file = self.log_file_path()
        self.write_log("Processing started")

        self.write_log("Preparing inputs")
        keywords = {}
        LOGGER.debug("received inputs: " + ", ".join(request.inputs.keys()))
        for name, input_queue in request.inputs.items():
            LOGGER.debug(input_queue)
            values = []

            for input_var in input_queue:
                if isinstance(input_var, ComplexInput):
                    ds = self.try_opendap(input_var)

                    if name in ds.data_vars:
                        value = ds.data_vars[name]
                    elif "variable" in request.inputs:
                        value = ds.data_vars[request.inputs["variable"][0].data]
                    else:
                        for key, val in ds.data_vars.items():
                            value = val
                            break
                        else:
                            raise ValueError(f"Input not understood: `{input_var}`.")

                elif isinstance(input, LiteralInput):
                    LOGGER.debug(input_var.data)
                    value = input_var.data

                else:
                    raise ValueError(f"Input not understood: `{input}`.")

                values.append(value)

            keywords[name] = values.pop() if len(values) == 1 else values
            keywords.pop("variable", None)

        self.write_log("Running computation")
        LOGGER.debug(keywords)
        out = self.xci(**keywords)
        out_fn = os.path.join(self.workdir, f"out_{self.identifier}.nc")

        self.write_log("Writing the output netcdf")
        out.to_netcdf(out_fn)
        response.outputs["output"].file = out_fn

        self.write_log("Processing finished successfully")
        return response


def chunk_dataset(ds, max_size=1000000):
    """
    Ensure the chunked size of a xarray.Dataset is below a certain size.

    Cycle through the dimensions, divide the chunk size by 2 until criteria is met.
    """
    chunks = dict(ds.sizes)

    def chunk_size():
        return reduce(mul, chunks.values())

    for dim in cycle(chunks):
        if chunk_size() < max_size:
            break
        chunks[dim] = max(chunks[dim] // 2, 1)

    return chunks


def make_freq(name, default="YS", allowed=("YS", "MS", "QS-DEC", "AS-JUL")):
    return LiteralInput(
        name,
        "Frequency",
        abstract="Resampling frequency",
        data_type="string",
        min_occurs=0,
        max_occurs=1,
        default=default,
        allowed_values=allowed,
    )


def make_thresh(name, default, abstract=""):
    return LiteralInput(
        name,
        "Threshold value, including units.",
        abstract=abstract,
        data_type="string",
        min_occurs=0,
        max_occurs=1,
        default=default,
    )


def make_window(name, default, abstract=""):
    return LiteralInput(
        name,
        "Window size",
        abstract=abstract,
        data_type="integer",
        min_occurs=0,
        max_occurs=1,
        default=default,
    )


def make_nc_input(name):
    return ComplexInput(
        name,
        "Resource",
        abstract="NetCDF Files or archive (tar/zip) containing netCDF files.",
        metadata=[Metadata("Info")],
        min_occurs=1,
        max_occurs=1000,
        supported_formats=[FORMATS.NETCDF],
    )


def make_mode(name, abstract=""):
    return LiteralInput(
        name,
        "Mode",
        abstract=abstract,
        data_type="string",
        allowed_values=["min", "max"],
        min_occurs=1,
        max_occurs=1,
    )


def make_op(name, abstract=""):
    return LiteralInput(
        name,
        "Operation name",
        abstract=abstract,
        data_type="string",
        allowed_values=[
            "min",
            "max",
            "mean",
            "std",
            "var",
            "count",
            "sum",
            "argmax",
            "argmin",
        ],
        min_occurs=1,
        max_occurs=1,
    )


def make_t(name, abstract=""):
    return LiteralInput(
        name,
        "Return period",
        abstract=abstract,
        data_type="integer",
        min_occurs=1,
        max_occurs=100,
    )


# TODO: Set allowed values to scipy.dist
def make_dist(name, abstract=""):
    return LiteralInput(
        name,
        "Distribution",
        abstract=abstract,
        data_type="string",
        min_occurs=1,
        max_occurs=1,
    )


def make_indexer():
    return [
        LiteralInput(
            "season",
            "Season",
            abstract="Season selection specification.",
            allowed_values=["DJF", "MAM", "JJA", "SON"],
            data_type="string",
            min_occurs=0,
            max_occurs=1,
            default=None,
        ),
        LiteralInput(
            "month",
            "Month",
            abstract="Month selection specification",
            data_type="string",
            allowed_values=range(1, 13),
            min_occurs=0,
            max_occurs=12,
            default=None,
        ),
    ]


def make_variable():
    return LiteralInput(
        "variable",
        "Variable name",
        abstract="Name of variable to analyze in netCDF file.",
        data_type="string",
        min_occurs=0,
        max_occurs=1,
    )
