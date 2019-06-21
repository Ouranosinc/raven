#########################################################################
:FileType          rvp ASCII Raven 2.8.2
:WrittenBy         Juliane Mai & James Craig
:CreationDate      Sep 2018
#
# Emulation of HMETS simulation of Salmon River near Prince George
#------------------------------------------------------------------------
#

#-----------------------------------------------------------------
# Raven Properties file Template. Created by Raven v2.8.2
#-----------------------------------------------------------------

# tied parameters:
# (it is important for OSTRICH to find every parameter place holder somewhere in this file)
# (without this "par_x06" and "par_x10" wouldn't be detectable)
#    para_sum_x05_x06 = par_x05 + par_x06
#    para_sum_x09_x10 = par_x09 + par_x10

#-----------------------------------------------------------------
# Soil Classes
#-----------------------------------------------------------------
:SoilClasses
  :Attributes,
  :Units,
  TOPSOIL,
  PHREATIC,
:EndSoilClasses

#-----------------------------------------------------------------
# Land Use Classes
#-----------------------------------------------------------------
:LandUseClasses,
  :Attributes,        IMPERM,    FOREST_COV,
       :Units,          frac,          frac,
       FOREST,           0.0,           1.0,
:EndLandUseClasses

#-----------------------------------------------------------------
# Vegetation Classes
#-----------------------------------------------------------------
:VegetationClasses,
  :Attributes,        MAX_HT,       MAX_LAI, MAX_LEAF_COND,
       :Units,             m,          none,      mm_per_s,
       FOREST,             4,             5,             5,
:EndVegetationClasses

#-----------------------------------------------------------------
# Soil Profiles
#-----------------------------------------------------------------
:SoilProfiles
         LAKE, 0
         ROCK, 0
  DEFAULT_P, 2, TOPSOIL,    par_x20, PHREATIC,    par_x21,
# DEFAULT_P, 2, TOPSOIL, x(20)/1000, PHREATIC, x(21)/1000,
:EndSoilProfiles

#-----------------------------------------------------------------
# Global Parameters
#-----------------------------------------------------------------
:GlobalParameter         SNOW_SWI_MIN par_x09           # x(9)
:GlobalParameter         SNOW_SWI_MAX par_sum_x09_x10   # x(9)+x(10)
:GlobalParameter     SWI_REDUCT_COEFF par_x11           # x(11)
:GlobalParameter             SNOW_SWI 0.05   #not sure why/if needed...

#-----------------------------------------------------------------
# Soil Parameters
#-----------------------------------------------------------------
:SoilParameterList
  :Parameters,        POROSITY,      PERC_COEFF,  PET_CORRECTION, BASEFLOW_COEFF
       :Units,               -,             1/d,               -,            1/d
      TOPSOIL,             1.0,         par_x17,         par_x15,        par_x18
     PHREATIC,             1.0,             0.0,             0.0,        par_x19
 #    TOPSOIL,             1.0,           x(17),           x(15),          x(18)
 #   PHREATIC,             1.0,             0.0,             0.0,          x(19)
:EndSoilParameterList

#-----------------------------------------------------------------
# Land Use Parameters
#-----------------------------------------------------------------
:LandUseParameterList
  :Parameters, MIN_MELT_FACTOR, MAX_MELT_FACTOR,    DD_MELT_TEMP,  DD_AGGRADATION, REFREEZE_FACTOR,    REFREEZE_EXP, DD_REFREEZE_TEMP, HMETS_RUNOFF_COEFF,
       :Units,          mm/d/C,          mm/d/C,               C,            1/mm,          mm/d/C,               -,                C,                  -,
    [DEFAULT],         par_x05, par_sum_x05_x06,         par_x07,         par_x08,         par_x13,         par_x14,          par_x12,            par_x16,
#                         x(5),       x(5)+x(6),            x(7),            x(8),           x(13),           x(14),            x(12),              x(16),
:EndLandUseParameterList
:LandUseParameterList
  :Parameters,   GAMMA_SHAPE,     GAMMA_SCALE,    GAMMA_SHAPE2,    GAMMA_SCALE2,
       :Units,             -,               -,               -,               -,
    [DEFAULT],       par_x01,         par_x02,         par_x03,         par_x04,
    #                   x(1),            x(2),            x(3),            x(4),
:EndLandUseParameterList
#-----------------------------------------------------------------
# Vegetation Parameters
#-----------------------------------------------------------------
:VegetationParameterList
  :Parameters,  RAIN_ICEPT_PCT,  SNOW_ICEPT_PCT,
       :Units,               -,               -,
    [DEFAULT],             0.0,             0.0,
:EndVegetationParameterList
