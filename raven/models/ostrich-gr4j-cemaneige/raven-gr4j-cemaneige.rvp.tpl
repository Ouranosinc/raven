#########################################################################
:FileType          rvp ASCII Raven 2.8.2
:WrittenBy         Juliane Mai & James Craig
:CreationDate      Sep 2018
#
# Emulation of GR4J simulation of Salmon River near Prince George
#------------------------------------------------------------------------

# tied parameters:
# (it is important for OSTRICH to find every parameter place holder somewhere in this file)
# (without this "para_x6" wouldn't be detectable)
#    para_1_minus_x6 = par_x6 - 1.0 = par_1_minus_x6

# -Global snow parameters-------------------------------------
:RainSnowTransition 0 1.0
:AirSnowCoeff       par_1_minus_x6  # [1/d] = 1.0 - CEMANEIGE_X2 = 1.0 - x6
:AvgAnnualSnow      par_x5          # [mm]  =       CEMANEIGE_X1 =       x5

# -Orographic Corrections-------------------------------------
:PrecipitationLapseRate 0.0004
:AdiabaticLapseRate 0.0065

# - Soil classes ---------------------------------------------
:SoilClasses
  :Attributes
  :Units
   SOIL_PROD
   SOIL_ROUT
   SOIL_TEMP
   SOIL_GW
:EndSoilClasses
:SoilParameterList
 :Parameters, POROSITY ,  GR4J_X3, GR4J_X2
 :Units     ,     none ,       mm,    mm/d
   [DEFAULT],      1.0 ,   par_x3,  par_x2
:EndSoilParameterList

# ----Soil Profiles--------------------------------------------
#     name, #horizons, (soiltype, thickness) x #horizons
#     GR4J_X1 is thickness of first layer (SOIL_PROD), here 0.696
:SoilProfiles
  DEFAULT_P,      4, SOIL_PROD   , par_x1, SOIL_ROUT  ,   0.300, SOIL_TEMP  ,   1.000, SOIL_GW  ,   1.000,
:EndSoilProfiles

# ----Vegetation Classes---------------------------------------
:VegetationClasses
   :Attributes,       MAX_HT,       MAX_LAI,      MAX_LEAF_COND
        :Units,            m,          none,           mm_per_s
       VEG_ALL,           0.0,          0.0,                0.0
:EndVegetationClasses

# --Land Use Classes------------------------------------------
:LandUseClasses
  :Attributes, IMPERM, FOREST_COV
  :Units     ,   frac,       frac
       LU_ALL,    0.0,        0.0
:EndLandUseClasses
:LandUseParameterList
 :Parameters, GR4J_X4, MELT_FACTOR
 :Units     ,       d,      mm/d/C
   [DEFAULT],  par_x4,        7.73
:EndLandUseParameterList
