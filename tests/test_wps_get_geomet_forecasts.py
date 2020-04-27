import json

from pywps import Service
from pywps.tests import assert_response_success
from shapely.geometry import shape, MultiPolygon

from raven.processes import GetGeometDataProcess
from .common import client_for, TESTDATA, CFG_FILE, get_output, count_pixels
import datetime as dt

class TestGenericZonalStatsProcess:

    def test_simple(self):
        client = client_for(Service(processes=[GetGeometDataProcess(), ], cfgfiles=CFG_FILE))

        fields = [
            'use_rdps={use_rdps}',
            'use_reps={use_reps}',
            'use_gdps={use_gdps}',
            'use_geps={use_geps}',
            'combine_reg_glob={combine_reg_glob}',
            'shape=file@xlink:href=file://{shape}',
            'forecast_date={forecast_date}',
        ]

        datainputs = ';'.join(fields).format(
            shape=TESTDATA['basin_test'],
            use_rdps=True,
            use_reps=True,
            use_gdps=True,
            use_geps=True,
            combine_reg_glob=True,
            forecast_date = dt.datetime.today()
        )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='geomet-forecasts', datainputs=datainputs)

        assert_response_success(resp)
        out = get_output(resp.xml)
