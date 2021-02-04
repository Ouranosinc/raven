"""
This module contains the WPS inputs and outputs that are reused across multiple WPS processes.


"""

from pywps import (
    FORMATS,
    ComplexInput,
    ComplexOutput,
    Format,
    LiteralInput,
    LiteralOutput,
)
from pywps.app.Common import Metadata
from ravenpy.models import GR4JCN, HBVEC, HMETS, MOHYSE, rv

from raven import config

# ---------------------------------------- #
# ---------------- Inputs ---------------- #
# ---------------------------------------- #


ts = ComplexInput(
    "ts",
    "Input time series files",
    abstract="Files (text or netCDF) storing"
    "daily liquid precipitation (pr), "
    "solid precipitation (prsn), "
    "minimum temperature (tasmin), "
    "maximum temperature (tasmax), "
    "potential evapotranspiration (evspsbl) and "
    "observed streamflow (qobs [m3/s]).",
    min_occurs=1,
    max_occurs=100,
    supported_formats=[FORMATS.NETCDF, FORMATS.DODS, FORMATS.TEXT, FORMATS.SHP],
)

conf = ComplexInput(
    "conf",
    "Raven/Ostrich configuration files",
    abstract="Model configuration files, including the primary input file (rvi), the parameter "
    "input file (rvp), the basin definition file (rvh), the time series input file "
    "(rvt), the initial conditions file (rvc). For Ostrich, include the Ostrich "
    "calibration config (txt) and templates (tpl).",
    min_occurs=5,
    max_occurs=5,
    supported_formats=[FORMATS.TEXT],
)

rvi = ComplexInput(
    "rvi",
    "Primary input file",
    abstract="The primary input file stores the model simulation options and numerical options.",
    min_occurs=1,
    max_occurs=1,
    supported_formats=[FORMATS.TEXT],
)

rvp = ComplexInput(
    "rvp",
    "Classed parameter input file",
    abstract="The classed parameter input file stores a database of soil, vegetation, river, "
    "aquifer, and land class pro-perties. Not all classes specified in the *.rvp file "
    "need to be included in the model.",
    min_occurs=1,
    max_occurs=1,
    supported_formats=[FORMATS.TEXT],
)

rvh = ComplexInput(
    "rvh",
    "HRU / Basin definition file",
    abstract="The HRU/basin definition file describes the topology of the basin network and the "
    "class membership of all constituent HRUs.",
    min_occurs=1,
    max_occurs=1,
    supported_formats=[FORMATS.TEXT],
)

rvt = ComplexInput(
    "rvt",
    "Time series input file",
    abstract="The time series input file is used to store time series of forcing functions ("
    "precipitation, temperature, etc.).",
    min_occurs=1,
    max_occurs=1,
    supported_formats=[FORMATS.TEXT],
)

rvc = ComplexInput(
    "rvc",
    "Initial conditions input file",
    abstract="The initial conditions input file is used to store the initial conditions for the "
    "model. By default, the initial conditions for all model state variables is zero, "
    "and there are no required commands in this file (it could even be completely "
    "empty).",
    min_occurs=0,
    supported_formats=[FORMATS.TEXT],
)

start_date = LiteralInput(
    "start_date",
    "Simulation start date (AAAA-MM-DD)",
    abstract="Start date of the simulation (AAAA-MM-DD). "
    "Defaults to the start of the forcing file. ",
    data_type="dateTime",
    default="0001-01-01 00:00:00",
    min_occurs=0,
    max_occurs=config.max_parallel_processes,
)

end_date = LiteralInput(
    "end_date",
    "Simulation end date (AAAA-MM-DD)",
    abstract="End date of the simulation (AAAA-MM-DD). "
    "Defaults to the end of the forcing file.",
    data_type="dateTime",
    default="0001-01-01 00:00:00",
    min_occurs=0,
    max_occurs=config.max_parallel_processes,
)

duration = LiteralInput(
    "duration",
    "Simulation duration (days)",
    abstract="Number of simulated days, defaults to the length of the input forcings.",
    data_type="nonNegativeInteger",
    default=0,
    min_occurs=0,
    max_occurs=config.max_parallel_processes,
)

run_name = LiteralInput(
    "run_name",
    "Simulation name",
    abstract="The name given to the simulation, for example <watershed>_<experiment>",
    data_type="string",
    default="raven-gr4j-cemaneige-sim",
    min_occurs=0,
    max_occurs=config.max_parallel_processes,
)

name = LiteralInput(
    "name",
    "Watershed name",
    abstract="The name of the watershed the model is run for.",
    data_type="string",
    default="watershed",
    min_occurs=0,
    max_occurs=config.max_parallel_processes,
)

area = LiteralInput(
    "area",
    "Watershed area (km2)",
    abstract="Watershed area (km2)",
    data_type="float",
    default=0.0,
    min_occurs=0,
    max_occurs=config.max_parallel_processes,
)

latitude = LiteralInput(
    "latitude",
    "Latitude",
    abstract="Watershed's centroid latitude",
    data_type="float",
    min_occurs=1,
    max_occurs=config.max_parallel_processes,
)

longitude = LiteralInput(
    "longitude",
    "Longitude",
    abstract="Watershed's centroid longitude",
    data_type="float",
    min_occurs=1,
    max_occurs=config.max_parallel_processes,
)

elevation = LiteralInput(
    "elevation",
    "Elevation (m)",
    abstract="Watershed's mean elevation (m)",
    data_type="float",
    min_occurs=1,
    max_occurs=config.max_parallel_processes,
)

model_name = LiteralInput(
    "model_name",
    "Hydrological model identifier",
    abstract="Hydrological model identifier: {HMETS, GR4JCN, MOHYSE}",
    data_type="string",
    allowed_values=("HMETS", "GR4JCN", "MOHYSE"),
    min_occurs=1,
    max_occurs=config.max_parallel_processes,
)

nc_index = LiteralInput(
    "nc_index",
    "NetCDF site coordinate index",
    abstract="The site index for a multi-basin netCDF file. This is ONLY necessary if the "
    "NetCDF variable is 2-dimensional (time, site).",
    data_type="integer",
    min_occurs=0,
    max_occurs=config.max_parallel_processes,
)

suppress_output = LiteralInput(
    "suppress_output",
    "Do not write hydrograph to disk",
    abstract="If True (default), hydrographs are not written to disk and thus not"
    "returned.",
    data_type="boolean",
    default=True,
)

rain_snow_fraction = LiteralInput(
    "rain_snow_fraction",
    "Rain snow partitioning",
    abstract="Algorithm used to partition rain and snow from the total precipitions",
    data_type="string",
    allowed_values=rv.rain_snow_fraction_options,
    min_occurs=0,
)

evaporation = LiteralInput(
    "evaporation",
    "Evaporation scheme",
    abstract="Algorithm used to compute potential evapotranspiration (PET).",
    data_type="string",
    allowed_values=rv.evaporation_options,
    min_occurs=0,
)

ow_evaporation = LiteralInput(
    "ow_evaporation",
    "Open-water evaporation scheme",
    abstract="Algorithm used to compute potential evapotranspiration (PET) over open "
    "water",
    data_type="string",
    allowed_values=rv.evaporation_options,
    min_occurs=0,
)

nc_spec = LiteralInput(
    "nc_spec",
    "NetCDF input file specifications",
    abstract="Configuration of individual netCDF input files, such as `linear_transform`"
    "and `time_shift`. Should be passed as a dictionary keyed by variable, e.g. `tas` "
    "json-serialized.",
    data_type="string",
    min_occurs=0,
    max_occurs=20,
)

forecast_model = LiteralInput(
    "forecast_model",
    "ECCC forecast model",
    abstract="The name of the forecast model run by Environment and Climate Change "
    "Canada.",
    data_type="string",
    allowed_values=("GEPS",),  # 'REPS', 'GDPS', 'RDPS'),
    default="GEPS",
    min_occurs=1,
)

hdate = LiteralInput(
    "hdate",
    "Hindcast start date (AAAA-MM-DD)",
    abstract="Start date of the hindcast (AAAA-MM-DD). "
    "Defaults to the start of the forcing file. ",
    data_type="dateTime",
    min_occurs=1,
    max_occurs=1,
)

hmets = LiteralInput(
    "hmets",
    "Comma separated list of HMETS parameters",
    abstract="Parameters: " + ", ".join(HMETS.params._fields),
    data_type="string",
    min_occurs=0,
)

gr4jcn = LiteralInput(
    "gr4jcn",
    "Comma separated list of GR4JCN parameters",
    abstract="Parameters: " + ", ".join(GR4JCN.params._fields),
    data_type="string",
    min_occurs=0,
)

mohyse = LiteralInput(
    "mohyse",
    "Comma separated list of MOHYSE parameters",
    abstract="Parameters: " + ", ".join(MOHYSE.params._fields),
    data_type="string",
    min_occurs=0,
)

hbvec = LiteralInput(
    "hbvec",
    "Comma separated list of HBV-EC parameters",
    abstract="Parameters: " + ", ".join(HBVEC.params._fields),
    data_type="string",
    min_occurs=0,
)

# --- GIS Inputs --- #

region_vector = ComplexInput(
    "region_vector",
    "Vector shape file of a region",
    abstract="An ESRI Shapefile, GML, JSON, GeoJSON, or single layer GeoPackage."
    " The ESRI Shapefile must be zipped and contain the .shp, .shx, and .dbf.",
    min_occurs=1,
    max_occurs=1,
    supported_formats=[
        FORMATS.GEOJSON,
        FORMATS.GML,
        FORMATS.JSON,
        FORMATS.SHP,
        FORMATS.ZIP,
    ],
)

shape = ComplexInput(
    "shape",
    "Vector shape of a region",
    abstract="An ESRI Shapefile, GML, JSON, GeoJSON, or single layer GeoPackage."
    " The ESRI Shapefile must be zipped and contain the .shp, .shx, and .dbf.",
    min_occurs=1,
    max_occurs=1,
    supported_formats=[
        FORMATS.GEOJSON,
        FORMATS.GML,
        FORMATS.JSON,
        FORMATS.SHP,
        FORMATS.ZIP,
    ],
)

land_use_raster = ComplexInput(
    "raster",
    "Gridded Land Use raster data set",
    abstract="The Land Use raster to be queried. Default is the CEC NALCMS 2010. Provided "
    "raster "
    "must use the UN FAO Land Cover Classification System (19 types).",
    metadata=[
        Metadata(
            "Commission for Environmental Cooperation North American Land Change Monitoring "
            "System",
            "http://www.cec.org/tools-and-resources/map-files/land-cover-2010-landsat-30m",
        ),
        Metadata(
            "Latifovic, R., Homer, C., Ressl, R., Pouliot, D., Hossain, S.N., Colditz, R.R.,"
            "Olthof, I., Giri, C., Victoria, A., (2012). North American land change "
            "monitoring system. In: Giri, C., (Ed), Remote Sensing of Land Use and Land "
            "Cover: Principles and Applications, CRC-Press, pp. 303-324"
        ),
    ],
    min_occurs=0,
    max_occurs=1,
    supported_formats=[FORMATS.GEOTIFF],
)

dem_raster = ComplexInput(
    "raster",
    "Gridded raster data set",
    abstract="The DEM to be queried. Defaults to the EarthEnv-DEM90 product.",
    metadata=[
        Metadata("EarthEnv-DEM90", "https://www.earthenv.org/DEM"),
        Metadata(
            "Robinson, Natalie, James Regetz, and Robert P. Guralnick (2014). "
            "EarthEnv-DEM90: A Nearly-Global, Void-Free, Multi-Scale Smoothed, 90m Digital "
            "Elevation Model from Fused ASTER and SRTM Data. ISPRS Journal of "
            "Photogrammetry and Remote Sensing 87: 57–67.",
            "https://doi.org/10.1016/j.isprsjprs.2013.11.002",
        ),
    ],
    min_occurs=0,
    max_occurs=1,
    supported_formats=[FORMATS.GEOTIFF],
)

simple_categories = LiteralInput(
    "simple_categories",
    "Use simplified land classification categories for hydrological "
    "modeling purposes.",
    data_type="boolean",
    default="false",
    min_occurs=0,
    max_occurs=1,
)

raster_band = LiteralInput(
    "band",
    "Raster band",
    data_type="integer",
    default=1,
    abstract="Band of raster examined to perform zonal statistics.",
    min_occurs=0,
    max_occurs=1,
)

select_all_touching = LiteralInput(
    "select_all_touching",
    "Additionally select boundary pixels that are touched by shape.",
    data_type="boolean",
    default="false",
    min_occurs=0,
    max_occurs=1,
)


# --- #

rv_config = ComplexOutput(
    "rv_config",
    "Raven/Ostrich configuration files",
    abstract="Model configuration files, including the primary input file (rvi), the parameter "
    "input file (rvp), the basin definition file (rvh), the time series input file "
    "(rvt), the initial conditions file (rvc). For Ostrich, include the Ostrich "
    "calibration config (txt) and templates (tpl).",
    supported_formats=[FORMATS.ZIP],
    as_reference=True,
)

hydrograph = ComplexOutput(
    "hydrograph",
    "Hydrograph time series (m3/s)",
    supported_formats=[
        FORMATS.NETCDF,
        Format("application/zip", extension=".zip", encoding="base64"),
    ],
    abstract="A netCDF file containing the outflow hydrographs (in m3/s) for all subbasins "
    "specified as `gauged` in the .rvh file. It reports period-ending time-"
    "averaged flows for the preceding time step, as is consistent with most "
    "measured stream gauge data (again, the initial flow conditions at the "
    "start of the first time step are included). If observed hydrographs are "
    "specified, they will be output adjacent to the corresponding modelled  "
    "hydrograph. ",
    as_reference=True,
)

ensemble = ComplexOutput(
    "ensemble",
    "Multiple hydrograph time series (m3/s)",
    supported_formats=[FORMATS.NETCDF],
    abstract="A netCDF file containing the outflow hydrographs (in m3/s) for the basin "
    "on which the regionalization method has been applied. The number of outflow "
    "hydrographs is equal to the number of donors (ndonors) passed to the method. "
    "The average of these hydrographs (either using equal or Inverse-Distance Weights) "
    'is the hydrograph generated in "hydrograph".',
    as_reference=True,
)

forecast = ComplexOutput(
    "forecast",
    "Multiple forecasted hydrograph time series (m3/s)",
    supported_formats=[FORMATS.NETCDF],
    abstract="A netCDF file containing the outflow hydrographs (in m3/s) for the basin "
    "on which the forecasting method has been applied. The number of members "
    "(hydrographs) is equal to the number of input weather forecast members "
    "passed to the method. ",
    as_reference=True,
)

storage = ComplexOutput(
    "storage",
    "Watershed storage time series (mm)",
    abstract="A netCDF file describing the total storage of water (in mm) in all water "
    "storage compartments for each time step of the simulation. Mass balance "
    "errors, cumulative input (precipitation), and output (channel losses) are "
    "also included. Note that the precipitation rates in this file are "
    "period-ending, i.e., this is the precipitation rate for the time step "
    "preceding the time stamp; all water storage variables represent "
    "instantaneous reports of the storage at the time stamp indicate.",
    supported_formats=[
        FORMATS.NETCDF,
        Format("application/zip", extension=".zip", encoding="base64"),
    ],
    as_reference=True,
)

solution = ComplexOutput(
    "solution",
    "solution.rvc file to restart another simulation with the conditions "
    "at the end of this simulation.",
    supported_formats=[
        FORMATS.TEXT,
        Format("application/zip", extension=".zip", encoding="base64"),
    ],
    as_reference=True,
)

diagnostics = ComplexOutput(
    "diagnostics",
    "Performance diagnostic values",
    abstract="Model diagnostic CSV file.",
    supported_formats=[
        FORMATS.TEXT,
        Format("application/zip", extension=".zip", encoding="base64"),
    ],
    as_reference=True,
)

features = ComplexOutput(
    "features",
    "DEM properties within the region defined by the vector provided.",
    abstract="Category pixel counts using either standard or simplified UNFAO categories",
    supported_formats=[FORMATS.GEOJSON],
)
statistics = ComplexOutput(
    "statistics",
    "DEM properties by feature",
    abstract="Land-use type pixel counts using either standard or simplified UNFAO categories.",
    supported_formats=[FORMATS.JSON],
)

calibparams = LiteralOutput(
    "calibparams",
    "Calibrated prameters",
    abstract="Comma separated list of parameters.",
    data_type="string",
)

# --- OSTRICH --- #

algorithm = LiteralInput(
    "algorithm",
    "OSTRICH Algorithm to use to calibrate model parameters",
    abstract="Optimization algorithm to implement for this calibration run",
    data_type="string",
    default="DDS",
    allowed_values=("DDS", "SCEUA"),
    min_occurs=0,
)

max_iterations = LiteralInput(
    "max_iterations",
    "Maximum number of model evaluations for the calibration run (budget)",
    abstract="Maximum number of times OSTRICH can call the hydrological model during the "
    "model parameter calibrationn",
    data_type="integer",
    default=50,
    allowed_values=list(range(25001)),
    min_occurs=0,
)

random_seed = LiteralInput(
    "random_seed",
    "Seed for random number generator",
    abstract="Set this value to obtain replicable results. Set to -1 to let it be random.",
    data_type="integer",
    default=-1,
    min_occurs=0,
)

calibration = ComplexOutput(
    "calibration",
    "Ostrich calibration output",
    abstract="Output file from Ostrich calibration run.",
    supported_formats=[FORMATS.TEXT],
    as_reference=True,
)

CalibrationResults = ComplexOutput(
    "CalibrationResults",
    "ObjectiveFunction and calibrated parameters computed by Ostrich",
    abstract="Objective Function value after calibration using user-selected "
    "function, as well as the calibrated parameter set",
    supported_formats=[FORMATS.TEXT],
    as_reference=True,
)

calibrated_params = ComplexOutput(
    "calibrated_params",
    "Calibrated parameters",
    abstract="Model parameters estimated by minimizing the objective function.",
    supported_formats=[FORMATS.TEXT],
    as_reference=False,
)

# TODO: Add configuration files to output
# config = ComplexOutput('config', 'Configuration files',
#                        abstract="Link to configuration files.",
#                        supported_formats=)
