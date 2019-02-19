#########################################################################                                  
:FileType          rvc ASCII Raven 2.8.2                                                                              
:WrittenBy         Juliane Mai & James Craig                                                                             
:CreationDate      Sep 2018
#
# Emulation of GR4J simulation of Salmon River near Prince George                                                             
#------------------------------------------------------------------------                                 
#
# tied parameters:
# (it is important for OSTRICH to find every parameter place holder somewhere in this file)
# (without this "para_x1" wouldn't be detectable)
#    para_half_x1 = para_x1 / 2. = par_x1 / 2. [m] = par_half_x1 [mm]

# initialize to 1/2 full
# GR4J_X1 * 1000. / 2.0
:UniformInitialConditions SOIL[0] par_half_x1

# Fixed because SOIL_ROUT layer thickness is fixed to be 0.3m
:UniformInitialConditions SOIL[1] 15.0

:HRUStateVariableTable (formerly :InitialConditionsTable)
   :Attributes SOIL[0]       SOIL[1]
   :Units      mm            mm
   1           par_half_x1   15.0
:EndHRUStateVariableTable