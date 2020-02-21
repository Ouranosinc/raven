#########################################################################                                  
:FileType          rvh ASCII Raven 2.8.2                                                                              
:WrittenBy         Juliane Mai & James Craig                                                                             
:CreationDate      Nov 2018
#
# Emulation of HBV-EC simulation of Salmon River near Prince George                                                             
#------------------------------------------------------------------------   

#                                                                                                               
:SubBasins                                                                                                              
        :Attributes     NAME    DOWNSTREAM_ID   PROFILE   REACH_LENGTH    GAUGED                                                          
        :Units          none    none            none      km              none                                                                                                    
        1,               hbv,   -1,             NONE,     _AUTO,          1
:EndSubBasins                                                                                                                           
                                                                                                                
:HRUs                                                                                                           
        :Attributes     AREA    ELEVATION  LATITUDE    LONGITUDE   BASIN_ID  LAND_USE_CLASS  VEG_CLASS   SOIL_PROFILE  AQUIFER_PROFILE   TERRAIN_CLASS   SLOPE   ASPECT  
        :Units           km2            m       deg          deg       none            none       none           none             none            none   ratio      deg     
             1,       4250.6,       843.0,  54.4848,   -123.3659,         1,         LU_ALL,   VEG_ALL,     DEFAULT_P,          [NONE],         [NONE], [NONE],  [NONE]
:EndHRUs

:SubBasinProperties
#                       HBV_PARA_11, DERIVED FROM HBV_PARA_11,
#                            MAXBAS,                 MAXBAS/2,
   :Parameters,           TIME_CONC,             TIME_TO_PEAK 
   :Units     ,                   d,                        d,    
             1,             par_x11,             par_half_x11,    
:EndSubBasinProperties