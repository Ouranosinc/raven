#########################################################################                                  
:FileType          rvh ASCII Raven 2.8.2                                                                              
:WrittenBy         Juliane Mai & James Craig                                                                             
:CreationDate      Sep 2018
#
# Emulation of MOHYSE simulation of Salmon River near Prince George                                                             
#------------------------------------------------------------------------

# tied parameters:
# (it is important for OSTRICH to find every parameter place holder somewhere in this file)
# (without this "par_x11" wouldn't be detectable)
#    para_rezi_x10 = 1.0 / par_x10 
#    para_x11      = par_x11

:SubBasins                                                                                                              
        :Attributes     NAME    DOWNSTREAM_ID   PROFILE   REACH_LENGTH    GAUGED                                                          
        :Units          none    none            none      km              none                                                                                                    
        1,            mohyse,   -1,             NONE,     _AUTO,          1
:EndSubBasins                                                                                                                           
                                                                                                                
:HRUs                                                                                                           
        :Attributes     AREA    ELEVATION  LATITUDE    LONGITUDE   BASIN_ID  LAND_USE_CLASS  VEG_CLASS   SOIL_PROFILE  AQUIFER_PROFILE   TERRAIN_CLASS   SLOPE   ASPECT  
        :Units           km2            m       deg          deg       none            none       none           none             none            none   ratio      deg     
             1,       4250.6,       843.0,  54.4848,   -123.3659,         1,         LU_ALL,   VEG_ALL,     DEFAULT_P,          [NONE],         [NONE], [NONE],  [NONE]
:EndHRUs

:SubBasinProperties
#          1.0 / MOHYSE_PARA_10,   MOHYSE_PARA_9
   :Parameters,     GAMMA_SCALE,     GAMMA_SHAPE,
   :Units,                  1/d,               -
              1,   par_rezi_x10,         par_x09
:EndSubBasinProperties

