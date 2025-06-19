from pathlib import Path

from pywps import get_ElementMakerForVersion
from pywps.app.basic import get_xpath_ns
from pywps.tests import WpsClient, WpsTestResponse

VERSION = "1.0.0"
WPS, OWS = get_ElementMakerForVersion(VERSION)
xpath_ns = get_xpath_ns(VERSION)

CFG_FILE = [
    str(Path(__file__).parent / "test.cfg"),
]


class WpsTestClient(WpsClient):
    def get(self, *args, **kwargs):
        query = "?"
        for key, value in kwargs.items():
            query += f"{key}={value}&"
        return super().get(query)


def count_pixels(stats, numeric_categories=False):
    category_counts = 0
    for key, val in stats.items():
        if numeric_categories:
            try:
                int(key)
            except ValueError:
                continue
        if key in ["count", "min", "max", "mean", "median", "sum", "nodata"]:
            continue
        category_counts += val
    return category_counts


def client_for(service):
    return WpsTestClient(service, WpsTestResponse)


def get_output(doc):
    """Copied from pywps/tests/test_execute.py.
    TODO: make this helper method public in pywps."""
    output = {}
    for output_el in xpath_ns(
        doc, "/wps:ExecuteResponse/wps:ProcessOutputs/wps:Output"
    ):
        [identifier_el] = xpath_ns(output_el, "./ows:Identifier")

        lit_el = xpath_ns(output_el, "./wps:Data/wps:LiteralData")
        if lit_el:
            output[identifier_el.text] = lit_el[0].text

        ref_el = xpath_ns(output_el, "./wps:Reference")
        if ref_el:
            output[identifier_el.text] = ref_el[0].attrib["href"]

        data_el = xpath_ns(output_el, "./wps:Data/wps:ComplexData")
        if data_el:
            output[identifier_el.text] = data_el[0].text

    return output
