#########################################################################
:FileType          rvt ASCII Raven 2.8.2
:WrittenBy         Juliane Mai & James Craig
:CreationDate      Sep 2018
#
# Emulation of GR4J simulation
#-----------------------------

:Gauge meteorological forcings
   :Latitude    {latitude}
   :Longitude   {longitude}
   :Elevation   {elevation}
   :RainCorrection         par_x20 # HBV_PAR_20 == RFCF
   :SnowCorrection         par_x21 # HBV_PAR_21 == SFCF

   {monthly_ave_evaporation}
   {monthly_ave_temperature}

   {pr}
   {rainfall}
   {prsn}
   {tasmin}
   {tasmax}
   {tas}
   {evspsbl}
:EndGauge

# observed streamflow
{water_volume_transport_in_river_channel}
