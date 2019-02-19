#########################################################################                                  
:FileType          rvp ASCII Raven 2.8.2                                                                              
:WrittenBy         Juliane Mai & James Craig                                                                             
:CreationDate      Nov 2018
#
# Emulation of HBV simulation of Salmon River near Prince George                                                             
#------------------------------------------------------------------------

# tied parameters:
# (it is important for OSTRICH to find every parameter place holder somewhere in this file)
# (without this "para_x05" and "para_x15" wouldn't be detectable)
#    para_1_+_x15       = 1.0 + par_x15 
#    para_half_x11      = 0.5 * par_x11

#------------------------------------------------------------------------
# Global parameters
#
#                             HBV_PARA_13=TCALT
:AdiabaticLapseRate                     par_x13
#                                   HBV_PARA_01, CONSTANT,
:RainSnowTransition                     par_x01       2.0
#                                   HBV_PARA_04,
:IrreducibleSnowSaturation              par_x04
#                             HBV_PARA_12=PCALT
:GlobalParameter PRECIP_LAPSE           par_x12

#---------------------------------------------------------
# Soil classes
:SoilClasses
 :Attributes, 
 :Units,      
   TOPSOIL,      1.0,    0.0,       0
   SLOW_RES,     1.0,    0.0,       0
   FAST_RES,     1.0,    0.0,       0  
:EndSoilClasses

:SoilParameterList
  :Parameters,                POROSITY,FIELD_CAPACITY,    SAT_WILT,    HBV_BETA, MAX_CAP_RISE_RATE,MAX_PERC_RATE,BASEFLOW_COEFF,            BASEFLOW_N
  :Units     ,                    none,          none,        none,        none,              mm/d,         mm/d,           1/d,                  none
  #                        HBV_PARA_05,   HBV_PARA_06, HBV_PARA_14, HBV_PARA_07,       HBV_PARA_16,     CONSTANT,      CONSTANT,              CONSTANT,
    [DEFAULT],                 par_x05,       par_x06,     par_x14,     par_x07,           par_x16,          0.0,           0.0,                   0.0
  #                                                       CONSTANT,                                  HBV_PARA_08,   HBV_PARA_09, 1+HBV_PARA_15=1+ALPHA,                 
     FAST_RES,                _DEFAULT,      _DEFAULT,         0.0,    _DEFAULT,          _DEFAULT,      par_x08,       par_x09,           par_1_+_x15
  #                                                       CONSTANT,                                                 HBV_PARA_10,              CONSTANT,
     SLOW_RES,                _DEFAULT,      _DEFAULT,         0.0,    _DEFAULT,          _DEFAULT,     _DEFAULT,       par_x10,                   1.0
:EndSoilParameterList

#---------------------------------------------------------
# Soil profiles
# name, layers, (soilClass, thickness) x layers
#
:SoilProfiles
#                        HBV_PARA_17,           CONSTANT,           CONSTANT,
   DEFAULT_P, 3, TOPSOIL,    par_x17, FAST_RES,    100.0, SLOW_RES,    100.0
:EndSoilProfiles

#---------------------------------------------------------
# Vegetation classes
#
:VegetationClasses
 :Attributes,   MAX_HT,  MAX_LAI, MAX_LEAF_COND
 :Units,             m,     none,      mm_per_s
   VEG_ALL,         25,      6.0,           5.3
:EndVegetationClasses

:VegetationParameterList
  :Parameters,  MAX_CAPACITY, MAX_SNOW_CAPACITY,  TFRAIN,  TFSNOW,
  :Units,                 mm,                mm,    frac,    frac,
  VEG_ALL,             10000,             10000,    0.88,    0.88,
:EndVegetationParameterList

#---------------------------------------------------------
# LandUse classes
#
:LandUseClasses
 :Attributes,     IMPERM, FOREST_COV
 :Units,            frac,       frac
      LU_ALL,        0.0,          1
:EndLandUseClasses

:LandUseParameterList
  :Parameters,   MELT_FACTOR, MIN_MELT_FACTOR,   HBV_MELT_FOR_CORR, REFREEZE_FACTOR, HBV_MELT_ASP_CORR
  :Units     ,        mm/d/K,          mm/d/K,                none,          mm/d/K,              none
  #              HBV_PARA_02,        CONSTANT,         HBV_PARA_18,     HBV_PARA_03,          CONSTANT
    [DEFAULT],       par_x02,             2.2,             par_x18,         par_x03,              0.48
:EndLandUseParameterList

:LandUseParameterList
 :Parameters, HBV_MELT_GLACIER_CORR,   HBV_GLACIER_KMIN, GLAC_STORAGE_COEFF, HBV_GLACIER_AG
 :Units     ,                  none,                1/d,                1/d,           1/mm  
   #                       CONSTANT,           CONSTANT,        HBV_PARA_19,       CONSTANT,
   [DEFAULT],                  1.64,               0.05,            par_x19,           0.05
:EndLandUseParameterList