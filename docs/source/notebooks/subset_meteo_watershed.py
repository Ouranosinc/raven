# %matplotlib inline

from birdy import WPSClient
import json
import xarray as xr
import matplotlib.pyplot


url_raven = "http://localhost:9099/wps"
url_flyingpigeon = 'http://localhost:8093/wps'

geoserver = 'https://boreas.ouranos.ca/geoserver/wfs'

nc1 = '/example/tas_Amon_CanESM2_rcp85_r1i1p1_200601-210012.nc'
nc2 = '<netCDF served via THREDDS>'

raven = WPSClient(url_raven, progress=True)  # processes='base_flow_index')
fp = WPSClient(url_flyingpigeon, progress=True)

print(raven, fp)

basin_process = raven.hydrosheds_select(location="-68.724444, 50.646667", aggregate_upstream=True)
feature, upstream_ids = basin_process.get(asobj=True)
featureids = upstream_ids[0:3]
# featureids should be 94396, 94436, 94437
print(featureids)

with open('testfile.json', 'w') as f:
    json.dump(feature, f)

x = fp.subset_wfs_polygon(resource=nc1, typename='public:USGS_HydroBASINS_lake_na_lev12',
                          geoserver=geoserver, featureids=featureids, mosaic='False')
ncfile, meta = x.get(asobj=True)
ds = xr.open_dataset(ncfile)
print(ds.vars, ds.attrs)
