Changes
=======

0.14.1
------
* Update to RavenPy 0.7.3 (pin xclim version 0.28.1)

0.14
----

* Update to RavenPy 0.7.2
* Use new OWSlib WFS topological filters
* More informative install documentation
* Upgrade to PyWPS 4.4.5

0.13
----

* Update RavenPy to 0.5.1
* Remove the ``name`` (watershed name) from the WPS interface for Raven processes
* Add ``random_numbers`` WPS param to pass optional ``OstRandomNumbers.txt`` file to Ostrich processes
* Add error handlers for regionalisation and climatology processes

0.12.1
------

* Fix bug where the name of configuration files was used, while the client transmission of data does not carry the file name.
* Update notebooks
* Move draft notebooks to sandbox

0.12.0
------

* Update RavenPy to 0.4.2
* Migrate utilities to RavenPy
* Add notebook for advanced forecasting
* Add notebook for probabilistic flood assessment
* Option to skip slow tests
* Add climpred verification WPS service
* Pre-commit hooks
* Install from conda Raven and Ostrich libraries
* Support passing HRUs
* Use scale/offset instead of linear_transform
* Enable GitHub CI
* Fix broken notebooks
* Improve error reporting by including stack trace in error messages.


0.11.x
------

* Add processes to run hydrological simulations on ECCC GEPS forecasts/hindcasts
* Add process to create forecast graphic
* Add first basic data assimilation utilities
* Factor out extra project RavenPy (at version 0.2.2), using Raven 3.0.1
* Upgrade to xclim +0.23.0
* Upgrade to xarray +0.16.2
* Add configuration options: ``deaccumulate``
* Clean notebooks
* Pin RavenPy to 0.3.0
* Pin owslib to 0.21
* Fix RavenC binaries installation for deployment
* Move some tests to RavenPy
* Regionalization data is now bundled with RavenPy
* Upgrade and pin PyWPS to 4.4.1
* Factor out most GIS functions to RavenPy (0.3.0)
* Add ``nalcms-zonal-stats-raster`` process using ``pymetalink``
* Simplify documentation build environment.


0.10.x (2020-03-09) Oxford
--------------------------

* ``suppress_ouput`` also triggers ``:DontWriteWatershedStorage``
* Added support for ERA5 (hourly), NRCan and CANOPEX datasets
* Support linear transforms (unit changes)
* Calibration now uses :SuppressOutput by default
* Added options for rain_snow_fraction, evaporation and ow_evaporation
* Updated Raven version to 295
* Support passing shapes as zip files


0.9.x (2019-11-11)
------------------

* Return configuration files used to run model in a zip archive


0.8.x (2019-10-22)
------------------
* Added more documentation for users
* Fixed reprojection errors in GIS utilities
* Specified HydroBASINS in lieu of HydroSHEDS in processes
* Optimized memory usage in ReadTheDocs builds when using Sphinx autodoc by employing mock
* Cleaner GeoJSON outputs for many subsetting processes
* Employed ipyleaflets for notebook-based web-maps
* Run py.test on notebooks from local or remote server


0.7.x (2019-06-25)
------------------

* Regionalization database
* Graphics for frequency analysis
* Many new notebook tutorials
* Bug fixes


0.6.x (2019-06-05)
------------------

* Regionalization process allowing the estimation of parameters of ungauged watersheds
* Added time series analysis processes, including frequential analysis
* Added processes creating graphics
* GIS processes now use GeoServer capabilities
* Docker configuration


0.5.0 (2019-04-12)
------------------

* Added watershed geospatial analysis processes
  - Hydroshed basin selection (with upstream contributors)
  - Watershed properties
  - DEM property analysis
  - Land-use property analysis
* Added multi-parameter parallel simulations
* Added multi-model parallel simulations
* Added multi-bassin parallel simulations


0.4.0 (2019-03-12)
------------------

* Added model calibration processes using Ostrich
* Added support for launching a singularity image
* Added library functions for model regionalization


0.3.0 (2019-01-24)
------------------

* Adds process for MOHYSE emulator
* Adds process for HBV-EC emulator


0.2.0 (2018-11-29) Washington
-----------------------------

* Provides generic RAVEN framework configuration
* Process for GR4J-Cemaneige emulator
* Process for HMETS emulator
