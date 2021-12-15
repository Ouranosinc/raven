=========
Notebooks
=========

These notebooks demonstrate a few features of the PAVICS-Hydro service.

If you're unfamiliar with notebooks, note that typing `TAB` after an object will display a drop-down menu of the object's attributes and methods, and that you need to either press the "run" button in the top of the JupyterLab window on the PAVICS server, or "hit `CTRL-Enter` to run a *cell*. You can also type `?` after a function or method to display the corresponding help message.


Getting started - Tutorial
==========================

.. toctree::
   :maxdepth: 1

   01_Introduction_to_RavenPy
   02_Emulating_hydrological_models
   03_Extract_geographical_watershed_properties
   04_Delineating_watersheds
   05_Extracting_external_data
   06_Raven_calibration
   07_Bias_correction_of_CMIP6_data
   08_Making_and_using_hotstart_files
   09_Hydrological_impacts_of_climate_change
   11_Climatological_ESP_forecasting
   12_Performing_hindcasting_experiments

More complex workflows making use of the RavenWPS server on the PAVICS platform
===============================================================================


General modelling
-----------------

.. toctree::
   :maxdepth: 1

   Running_hydrological_models_on_a_remote_server
   Running_HMETS_with_CANOPEX_dataset
   Model_calibration
   time_series_analysis
   Perform_Regionalization
   Region_selection


Forecasting
-----------

.. toctree::
   :maxdepth: 1

   Climatological_ESP_forecasting
   Hydrological_realtime_forecasting
   Hydrological_hindcasting
   Assess_probabilistic_flood_risk
   Comparing_hindcasts_and_ESP_forecasts


Parallel simulations
--------------------

.. toctree::
   :maxdepth: 1

   Multimodel_simulations
   Parameter_ensemble_simulations
   Multiple_watersheds_simulation

Advanced workflows
------------------

.. toctree::
   :maxdepth: 1

   Subset_climate_data_over_watershed
   Full_process_example_1

.. _OWSLib: https://geopython.github.io/OWSLib/
.. _Birdy: https://birdy.readthedocs.io
