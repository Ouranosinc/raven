import pytest
import datetime as dt

from pywps import Service
from pywps.tests import assert_response_success

from . common import client_for, TESTDATA, CFG_FILE
from raven.processes import RavenGR4JCemaNeigeProcess


class TestRavenGR4JCemaNeigeProcess:

    def test_simple(self):
        client = client_for(Service(processes=[RavenGR4JCemaNeigeProcess(), ], cfgfiles=CFG_FILE))

        datainputs = "nc=files@xlink:href=file://{fn};" \
                     "params={params};" \
                     "start_date={start_date};" \
                     "end_date={end_date};" \
                     "init={init};" \
                     "name={name};" \
                     "area={area};" \
                     "latitude={latitude};" \
                     "longitude={longitude};" \
                     "elevation={elevation};" \
            .format(fn=TESTDATA['raven-gr4j-cemaneige-nc-ts'],
                    params='0.696, 0.7, 19.7, 2.09, 123.3, 0.75',
                    start_date=dt.datetime(2000, 1, 1),
                    end_date=dt.datetime(2002, 1, 1),
                    init='4,5',
                    name='Saumon',
                    run_name='test',
                    area='4250.6',
                    elevation='843.0',
                    latitude=54.4848,
                    longitude=-123.3659,
                    )

        resp = client.get(
            service='WPS', request='Execute', version='1.0.0', identifier='raven-gr4j-cemaneige',
            datainputs=datainputs)

        assert_response_success(resp)
