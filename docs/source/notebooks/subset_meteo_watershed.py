# %matplotlib inline

from birdy import WPSClient

url_raven = "http://localhost:9099/wps"
url_flyingpigeon = 'http://localhost:8093/wps'

geoserver = 'https://boreas.ouranos.ca/geoserver/wfs'

nc1 = '<locally-stored NetCDF>'
nc2 = '<netCDF served via THREDDS>'

raven = WPSClient(url_raven, progress=True)  # processes='base_flow_index')
fp = WPSClient(url_flyingpigeon, progress=True)

print(raven, fp)

basin_process = raven.hydrosheds_select(location="-68.724444, 50.646667", aggregate_upstream=True)
feature, upstream_ids = basin_process.get(asobj=True)
featureids = upstream_ids[0:3]
# featureids = [94396, 94436, 94437]

print(featureids)

x = fp.subset_wfs_polygon(resource=nc1, typename='public:USGS_HydroBASINS_lake_na_lev12',
                          geoserver=geoserver, featureids=featureids, mosaic='False')


