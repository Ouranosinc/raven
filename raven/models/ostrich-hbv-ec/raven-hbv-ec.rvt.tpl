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

   #                       Jan        Feb        Mar        Apr        May        Jun        Jul        Aug        Sep        Oct        Nov        Dec
   :MonthlyAveEvaporation, {mae.x01}, {mae.x02}, {mae.x03}, {mae.x04}, {mae.x05}, {mae.x06}, {mae.x07}, {mae.x08}, {mae.x09}, {mae.x10}, {mae.x11}, {mae.x12}
   :MonthlyAveTemperature, {mat.x01}, {mat.x02}, {mat.x03}, {mat.x04}, {mat.x05}, {mat.x06}, {mat.x07}, {mat.x08}, {mat.x09}, {mat.x10}, {mat.x11}, {mat.x12}

   :Data RAINFALL mm/d
      :ReadFromNetCDF
         :FileNameNC      {pr}
         :VarNameNC       {pr_var}
         :DimNamesNC      {nc_dimensions}
         :StationIdx      {nc_index}
         {time_shift}
         :LinearTransform {linear_transform}
      :EndReadFromNetCDF
   :EndData
   :Data SNOWFALL mm/d
      :ReadFromNetCDF
         :FileNameNC      {prsn}
         :VarNameNC       {prsn_var}
         :DimNamesNC      {nc_dimensions}
         :StationIdx      {nc_index}
         {time_shift}
         :LinearTransform {linear_transform}
      :EndReadFromNetCDF
   :EndData
   :Data TEMP_MIN deg_C
      :ReadFromNetCDF
         :FileNameNC      {tasmin}
         :VarNameNC       {tasmin_var}
         :DimNamesNC      {nc_dimensions}
         :StationIdx      {nc_index}
         {time_shift}
         :LinearTransform {linear_transform}
      :EndReadFromNetCDF
   :EndData
   :Data TEMP_MAX deg_C
      :ReadFromNetCDF
         :FileNameNC      {tasmax}
         :VarNameNC       {tasmax_var}
         :DimNamesNC      {nc_dimensions}
         :StationIdx      {nc_index}
         {time_shift}
         :LinearTransform {linear_transform}
      :EndReadFromNetCDF
   :EndData
   :Data PET mm/d
      :ReadFromNetCDF
         :FileNameNC      {evspsbl}
         :VarNameNC       {evspsbl_var}
         :DimNamesNC      {nc_dimensions}
         :StationIdx      {nc_index}
         {time_shift}
         :LinearTransform {linear_transform}
      :EndReadFromNetCDF
   :EndData
:EndGauge

# observed streamflow
:ObservationData HYDROGRAPH 1 m3/s
    :ReadFromNetCDF
         :FileNameNC      {water_volume_transport_in_river_channel}
         :VarNameNC       {water_volume_transport_in_river_channel_var}
         :DimNamesNC      {nc_dimensions}
         :StationIdx      {nc_index}
         {time_shift}
         :LinearTransform {linear_transform}
    :EndReadFromNetCDF
:EndObservationData

