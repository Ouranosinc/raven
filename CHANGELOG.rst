Changelog
=========

v0.19.0 (2025-07-16)
--------------------

Major Changes
^^^^^^^^^^^^^
* The `raven` library is now available for installation on PyPI (under the name `birdhouse-raven`) (#536, #569)
* Documentation has been overhauled so that WPS processes are listed in a more organized manner (#535)
* Reimplemented the `cruft` configuration and ported recent changes to ensure that the project is up-to-date (#536)
* Removed `pysheds` and watershed delineation process due to a licensing issue (#495)
* Updated several dependencies to ensure that processes emit fewer warnings and notebooks render properly (#535)

Internal Changes
^^^^^^^^^^^^^^^^
* Allowed for larger int fields when writing to GeoJSON (#535)
* Ensured that JSON-serialized output does not include `numpy` complex types (#535)
* Added a GitHub Workflow to test the Dockerfile recipe configuration for RavenWPS (#481)
* Cleaned up the Dockerfile recipe configuration for `raven`. Now using `gunicorn` for service management (#481)
* Testing data fetching mechanism has been refactored to use `pooch` for better maintainability (#569)
* Tooling has been updated to use `ruff` and other newer tools for code quality and formatting (#569)
* `raven` now uses Trusted Publisher for TestPyPI and PyPI releases (#569)

v0.18.2 (2023-07-06)
--------------------

* Removed pin on `owslib` below v0.29 and pin on `fiona` below v2.0 (#480)
* Added a GitHub Workflow to test against macOS builds (#480)
* Adapted zonal statistics processes to support the latest `fiona` and `zonalstats` API changes (#480)

v0.18.1 (2023-05-23)
--------------------

* Removed obsolete components related to HPC interfaces (#477)
* Added Python3.11 to supported versions with a build for Python3.11 in CI (#477)
* Adjusted ReadMe to reflect recent significant changes (#477)
* Updated deprecated GitHub Actions (#477)

v0.18.0 (2023-05-23)
--------------------

Major Changes
^^^^^^^^^^^^^
* `singularity` components have been removed from `raven` (#470)
* Removed Raven WPS capabilities for hydrological modelling, graphing and forecasting (moved to RavenPy) (#464)
* Removed notebooks and migrated to Ravenpy. Adapted them to the new Ravenpy configuration (#464)
* Removed all tests related to Raven WPS modelling (#464)
* Raise error message if shape area for NALCMS zonal stats raster is above 100,000 km2 (#473)

v0.17.1 (2023-04-04)
--------------------

Internal Changes
^^^^^^^^^^^^^^^^
* Dockerfile configuration now uses Python3.8 and `condaforge/mambaforge` base image (#466)
* `pandas` is temporarily pinned below v2.0 (#466)

v0.17.0 (2023-02-28)
--------------------

Major Changes
^^^^^^^^^^^^^
* Updated testing ensemble to use `pytest-xdist` (#448)
* Updated `RavenPy` to v0.11.0, `raven-hydro` to v3.6, and `fiona` to v1.9 (#461)
* Modified several geospatial processes to adapt to new APIs (#461)
* Datetime signatures for some models used in notebooks have been adjusted/fixed (#453)

Internal Changes
^^^^^^^^^^^^^^^^
* Makefile updates to better perform notebook refresh actions (#459)
* Pre-commit style updates (#446, #447, #449, #461)
* Use `provision-with-micromamba` GitHub Action in CI workflows (#462)

v0.16.0 (2022-11-01)
--------------------

Major Changes
^^^^^^^^^^^^^
* Added data assimilation workbook (#421)
* Overhaul of all existing notebooks within documentation (#424)
* Added notebooks for case-study paper (#435)
* Update to RavenPy 0.8.1 (#439)
* Dropped support for Python3.7 (#439)

Internal Changes
^^^^^^^^^^^^^^^^
* Added pre-commit.ci to workflows and updated black formatting (#428 and #429)
* Adjust documentation to remove sphinx-autoapi artefact files and set ReadTheDocs to fail_on_warning (#439)
* Set pre-commit to run new correction and verification hooks (#439):
    - pyupgrade: Ensure that coding style uses Python3.8+ conventions
    - pygrep: Checks for bare `noqa` comments and malformed code blocks in documentation
    - nbqa: Black, Isort, PyUpgrade now runs over notebooks
    - check-manifest: Ensure relevant modules and data are explicitly installed
    - black + blackdoc + yamllint: Clean up code, code examples within documentation and reformat yaml files for readability
    - check-jsonschema: Verify that GitHub and ReadTheDocs workflows are valid
* Added a Zenodo/DOI configuration

v0.15.1 (2022-01-14)
--------------------

* Modified handling for GDAL to better support conda-build configuration
* Update to RavenPy 0.7.8
* Upgrade to PyWPS 4.5.1

v0.15.0 (2021-12-22)
--------------------

* Update to RavenPy 0.7.7
* Update required Python consistently to v3.7+
* Set development status to Beta.
* Replace pip-installed packages with conda-forge equivalents.

v0.14.2 (2021-09-03)
--------------------

* Update to RavenPy 0.7.4 (pin climpred below version 2.1.6)
* Fixed a process-breaking bug in `wps_hydrobasins_shape_selection`

v0.14.1 (2021-08-31)
--------------------

* Update to RavenPy 0.7.3 (pin xclim version 0.28.1)

v0.14.0 (2021-08-30)
--------------------

* Update to RavenPy 0.7.2
* Use new OWSlib WFS topological filters
* More informative install documentation
* Upgrade to PyWPS 4.4.5

v0.13.0 (2021-05-14)
--------------------

* Update RavenPy to 0.5.1
* Remove the ``name`` (watershed name) from the WPS interface for Raven processes
* Add ``random_numbers`` WPS param to pass optional ``OstRandomNumbers.txt`` file to Ostrich processes
* Add error handlers for regionalisation and climatology processes

v0.12.1 (2021-04-16)
--------------------

* Fix bug where the name of configuration files was used, while the client transmission of data does not carry the file name.
* Update notebooks
* Move draft notebooks to sandbox

v0.12.0 (2021-04-14)
--------------------

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

v0.11.x (2021-02-01)
--------------------

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


v0.10.x (2020-03-09) Oxford
---------------------------

* ``suppress_ouput`` also triggers ``:DontWriteWatershedStorage``
* Added support for ERA5 (hourly), NRCan and CANOPEX datasets
* Support linear transforms (unit changes)
* Calibration now uses :SuppressOutput by default
* Added options for rain_snow_fraction, evaporation and ow_evaporation
* Updated Raven version to 295
* Support passing shapes as zip files

v0.9.x (2019-11-11)
-------------------

* Return configuration files used to run model in a zip archive

v0.8.x (2019-10-22)
--------------------
* Added more documentation for users
* Fixed reprojection errors in GIS utilities
* Specified HydroBASINS in lieu of HydroSHEDS in processes
* Optimized memory usage in ReadTheDocs builds when using Sphinx autodoc by employing mock
* Cleaner GeoJSON outputs for many subsetting processes
* Employed ipyleaflets for notebook-based web-maps
* Run py.test on notebooks from local or remote server

v0.7.x (2019-06-25)
-------------------

* Regionalization database
* Graphics for frequency analysis
* Many new notebook tutorials
* Bug fixes

v0.6.x (2019-06-05)
-------------------

* Regionalization process allowing the estimation of parameters of ungauged watersheds
* Added time series analysis processes, including frequential analysis
* Added processes creating graphics
* GIS processes now use GeoServer capabilities
* Docker configuration

v0.5.0 (2019-04-12)
-------------------

* Added watershed geospatial analysis processes
  - Hydroshed basin selection (with upstream contributors)
  - Watershed properties
  - DEM property analysis
  - Land-use property analysis
* Added multi-parameter parallel simulations
* Added multi-model parallel simulations
* Added multi-bassin parallel simulations

v0.4.0 (2019-03-12)
-------------------

* Added model calibration processes using Ostrich
* Added support for launching a singularity image
* Added library functions for model regionalization

v0.3.0 (2019-01-24)
-------------------

* Adds process for MOHYSE emulator
* Adds process for HBV-EC emulator

v0.2.0 (2018-11-29) Washington
------------------------------

* Provides generic RAVEN framework configuration
* Process for GR4J-Cemaneige emulator
* Process for HMETS emulator
