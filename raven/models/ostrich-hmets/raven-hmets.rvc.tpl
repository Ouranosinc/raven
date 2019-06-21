#########################################################################
:FileType          rvc ASCII Raven 2.8.2
:WrittenBy         Juliane Mai & James Craig
:CreationDate      Sep 2018
#
# Emulation of HMETS simulation of Salmon River near Prince George
#------------------------------------------------------------------------
#

# tied parameters:
# (it is important for OSTRICH to find every parameter place holder somewhere in this file)
# (without this "par_x20" and "par_x21" wouldn't be detectable)
#    para_half_x20 = para_x20 / 2. = par_x20 / 2. [m] = par_half_x20 [mm]
#    para_half_x21 = para_x21 / 2. = par_x21 / 2. [m] = par_half_x21 [mm]


# initialize to 1/2 full
:UniformInitialConditions SOIL[0] par_half_x20
:UniformInitialConditions SOIL[1] par_half_x21

:HRUStateVariableTable (formerly :IntialConditionsTable)
   :Attributes SOIL[0] SOIL[1]
   :Units mm mm
   1 par_half_x20 par_half_x21
:EndHRUStateVariableTable
