#########################################################################                                  
:FileType          rvt ASCII Raven 2.8.2                                                                              
:WrittenBy         Juliane Mai & James Craig                                                                             
:CreationDate      Sep 2018
#
# Emulation of GR4J simulation of Salmon River near Prince George                                                             
#------------------------------------------------------------------------

:Gauge meteorological forcings
   :Latitude    54.4848
   :Longitude -123.3659
   :Elevation  843.0
   :RainCorrection         par_x20 # 1.141608   # HBV_PAR_20 == RFCF
   :SnowCorrection         par_x21 # 1.024278   # HBV_PAR_21 == SFCF
   #                       Jan    Feb    Mar    Apr    May    Jun    Jul    Aug    Sep    Oct    Nov    Dec
   :MonthlyAveEvaporation, 0.028, 0.089, 0.315, 1.003, 2.077, 2.959, 3.190, 2.548, 1.382, 0.520, 0.101, 0.023,
   :MonthlyAveTemperature, -11.4, -7.2,  -2.8,  2.8,   8.1,   12.2,  14.4,  13.5,  8.9,   3.4,   -4.0,  -9.3,
   :Data RAINFALL mm/d
      :ReadFromNetCDF
         :FileNameNC     Salmon-River-Near-Prince-George_meteo_daily.nc
         :VarNameNC      rain
         :DimNamesNC     time
      :EndReadFromNetCDF
   :EndData
   :Data SNOWFALL mm/d
      :ReadFromNetCDF
         :FileNameNC     Salmon-River-Near-Prince-George_meteo_daily.nc
         :VarNameNC      snow
         :DimNamesNC     time
      :EndReadFromNetCDF
   :EndData
   :Data TEMP_MIN deg_C
      :ReadFromNetCDF
         :FileNameNC     Salmon-River-Near-Prince-George_meteo_daily.nc
         :VarNameNC      tmin
         :DimNamesNC     time
      :EndReadFromNetCDF
   :EndData
   :Data TEMP_MAX deg_C
      :ReadFromNetCDF
         :FileNameNC     Salmon-River-Near-Prince-George_meteo_daily.nc
         :VarNameNC      tmax
         :DimNamesNC     time
      :EndReadFromNetCDF
   :EndData
   :Data PET deg_C
      :ReadFromNetCDF
         :FileNameNC     Salmon-River-Near-Prince-George_meteo_daily.nc
         :VarNameNC      pet
         :DimNamesNC     time
      :EndReadFromNetCDF
   :EndData
:EndGauge

# # observed streamflow

:ObservationData	HYDROGRAPH	1	m3/s
    :ReadFromNetCDF 
         :FileNameNC     Salmon-River-Near-Prince-George_meteo_daily.nc
         :VarNameNC      qobs
         :DimNamesNC     time
    :EndReadFromNetCDF
:EndObservationData


