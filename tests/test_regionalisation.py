from raven.utilities import regionalization as reg
import datetime as dt
from .common import TESTDATA


def test_regionalization():
    reg.regionalize('SP', 'GR4JCN',
                    start_date=dt.datetime(2000, 1, 1),
                    end_date=dt.datetime(2002, 1, 1),
                    name='Salmon',
                    run_name='test',
                    area='4250.6',
                    elevation='843.0',
                    latitude=40.4848,
                    longitude=-103.3659,
                    min_NSE=0.6,
                    size=2,
                    ts=TESTDATA['raven-hmets-nc-ts'])
