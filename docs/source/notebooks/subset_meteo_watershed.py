# %matplotlib inline

from birdy import WPSClient
import json
import xarray as xr
import requests
import matplotlib.pyplot


url_raven = "http://localhost:9099/wps"
url_flyingpigeon = "http://localhost:8093/wps"

geoserver = "https://boreas.ouranos.ca/geoserver/wfs"

nc1 = "/home/tjs/Desktop/tas/tas_Amon_CanESM2_rcp85_r1i1p1_200601-210012.nc"
nc2 = "<netCDF served via THREDDS>"

raven = WPSClient(url_raven, progress=True)  # processes='base_flow_index')
fp = WPSClient(url_flyingpigeon, progress=True)

print(raven, fp)

basin_process = raven.hydrosheds_select(
    location="-68.724444, 50.646667", aggregate_upstream=True
)

url = basin_process.get(asobj=False)
feature = requests.get(url[0]).content
upstream_ids = json.loads(requests.get(url[1]).content)
featureids = upstream_ids[0:3]

with open("testfile.geojson", "wb") as f:
    f.write(feature)

# featureids should be 94396, 94436, 94437

x = fp.subset_wfs_polygon(
    resource=nc1,
    typename="public:USGS_HydroBASINS_lake_na_lev12",
    geoserver=geoserver,
    featureids=featureids,
    mosaic="False",
)

file_url, meta = x.get(asobj=False)

response = requests.get(file_url).content
with open('sdsdsd.nc', 'wb') as f:
    f.write(response)

ds = xr.open_dataset('sdsdsd.nc')
print(ds.variables)
