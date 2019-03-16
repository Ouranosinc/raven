#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 09:30:38 2019

@author: ets
"""

"""

BIRDY TEST

"""
from BirdyTestdata import TESTDATA
from birdy import WPSClient
import datetime as dt
from urllib.request import urlretrieve
import xarray as xr
import numpy as np
 
rwps = WPSClient("http://localhost:9099/wps")

params = '9.5019, 0.2774, 6.3942, 0.6884, 1.2875, 5.4134, 2.3641, 0.0973, 0.0464, 0.1998, 0.0222, -1.0919, ' \
                 '2.6851, 0.3740, 1.0000, 0.4739, 0.0114, 0.0243, 0.0069, 310.7211, 916.1947'
ts=TESTDATA['raven-hmets-nc-ts']
start_date=dt.datetime(2000, 1, 1)
end_date=dt.datetime(2002, 1, 1)
area=4250.6
elevation=843.0
latitude=54.4848
longitude=-123.3659

resp=rwps.raven_hmets(ts=str(ts),params=params,start_date=start_date, end_date=end_date, area=area, longitude=longitude, latitude=latitude, elevation=elevation)

[hydrograph,storage,solution,diagnostics]=resp.get(asobj=False)
tmp_file, _ = urlretrieve(diagnostics)
tmp_content = open(tmp_file).readlines()
idx_diag = tmp_content[0].split(',').index("DIAG_NASH_SUTCLIFFE")
diag = np.float(tmp_content[1].split(',')[idx_diag])

# Hydrograph
tmp_file, _ = urlretrieve(hydrograph)
ds = xr.open_dataset(tmp_file)


