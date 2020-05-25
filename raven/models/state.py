"""
State variables
===============

* HRUStateVariables
* BasinStateVariables

Implemented using typing.NamedTuple

Use _replace to update individual values.

"""

from typing import NamedTuple

_hru_state_variables_string = 'SURFACE_WATER,ATMOSPHERE,ATMOS_PRECIP,PONDED_WATER,SOIL[0],SOIL[1],SOIL[2],SOIL[3],' \
                              'SNOW_TEMP,SNOW,SNOW_COVER,AET,CONVOLUTION[0],CONVOLUTION[1],CONV_STOR[0],CONV_STOR[1],CONV_STOR[2],CONV_STOR[3],CONV_STOR[4],CONV_STOR[5],CONV_STOR[6],CONV_STOR[7],CONV_STOR[8],CONV_STOR[9],CONV_STOR[10],CONV_STOR[11],CONV_STOR[12],CONV_STOR[13],CONV_STOR[14],CONV_STOR[15],CONV_STOR[16],CONV_STOR[17],CONV_STOR[18],CONV_STOR[19],CONV_STOR[20],CONV_STOR[21],CONV_STOR[22],CONV_STOR[23],CONV_STOR[24],CONV_STOR[25],CONV_STOR[26],CONV_STOR[27],CONV_STOR[28],CONV_STOR[29],CONV_STOR[30],CONV_STOR[31],CONV_STOR[32],CONV_STOR[33],CONV_STOR[34],CONV_STOR[35],CONV_STOR[36],CONV_STOR[37],CONV_STOR[38],CONV_STOR[39],CONV_STOR[40],CONV_STOR[41],CONV_STOR[42],CONV_STOR[43],CONV_STOR[44],CONV_STOR[45],CONV_STOR[46],CONV_STOR[47],CONV_STOR[48],CONV_STOR[49],CONV_STOR[50],CONV_STOR[51],CONV_STOR[52],CONV_STOR[53],CONV_STOR[54],CONV_STOR[55],CONV_STOR[56],CONV_STOR[57],CONV_STOR[58],CONV_STOR[59],CONV_STOR[60],CONV_STOR[61],CONV_STOR[62],CONV_STOR[63],CONV_STOR[64],CONV_STOR[65],CONV_STOR[66],CONV_STOR[67],CONV_STOR[68],CONV_STOR[69],CONV_STOR[70],CONV_STOR[71],CONV_STOR[72],CONV_STOR[73],CONV_STOR[74],CONV_STOR[75],CONV_STOR[76],CONV_STOR[77],CONV_STOR[78],CONV_STOR[79],CONV_STOR[80],CONV_STOR[81],CONV_STOR[82],CONV_STOR[83],CONV_STOR[84],CONV_STOR[85],CONV_STOR[86],CONV_STOR[87],CONV_STOR[88],CONV_STOR[89],CONV_STOR[90],CONV_STOR[91],CONV_STOR[92],CONV_STOR[93],CONV_STOR[94],CONV_STOR[95],CONV_STOR[96],CONV_STOR[97],CONV_STOR[98],CONV_STOR[99]'

#_names = [v.replace('[', '').replace(']', '').lower() for v in _hru_state_variables_string.split(',')]


class HRUStateVariables(NamedTuple):
    surface_water: float = 0
    atmosphere: float = 0
    atmos_precip: float = 0
    ponded_water: float = 0
    soil0: float = 0
    soil1: float = 0
    soil2: float = 0
    soil3: float = 0
    snow_temp: float = 0
    snow: float = 0
    snow_cover: float = 0
    aet: float = 0
    convolution0: float = 0
    convolution1: float = 0
    conv_stor0: float = 0
    conv_stor1: float = 0
    conv_stor2: float = 0
    conv_stor3: float = 0
    conv_stor4: float = 0
    conv_stor5: float = 0
    conv_stor6: float = 0
    conv_stor7: float = 0
    conv_stor8: float = 0
    conv_stor9: float = 0
    conv_stor10: float = 0
    conv_stor11: float = 0
    conv_stor12: float = 0
    conv_stor13: float = 0
    conv_stor14: float = 0
    conv_stor15: float = 0
    conv_stor16: float = 0
    conv_stor17: float = 0
    conv_stor18: float = 0
    conv_stor19: float = 0
    conv_stor20: float = 0
    conv_stor21: float = 0
    conv_stor22: float = 0
    conv_stor23: float = 0
    conv_stor24: float = 0
    conv_stor25: float = 0
    conv_stor26: float = 0
    conv_stor27: float = 0
    conv_stor28: float = 0
    conv_stor29: float = 0
    conv_stor30: float = 0
    conv_stor31: float = 0
    conv_stor32: float = 0
    conv_stor33: float = 0
    conv_stor34: float = 0
    conv_stor35: float = 0
    conv_stor36: float = 0
    conv_stor37: float = 0
    conv_stor38: float = 0
    conv_stor39: float = 0
    conv_stor40: float = 0
    conv_stor41: float = 0
    conv_stor42: float = 0
    conv_stor43: float = 0
    conv_stor44: float = 0
    conv_stor45: float = 0
    conv_stor46: float = 0
    conv_stor47: float = 0
    conv_stor48: float = 0
    conv_stor49: float = 0
    conv_stor50: float = 0
    conv_stor51: float = 0
    conv_stor52: float = 0
    conv_stor53: float = 0
    conv_stor54: float = 0
    conv_stor55: float = 0
    conv_stor56: float = 0
    conv_stor57: float = 0
    conv_stor58: float = 0
    conv_stor59: float = 0
    conv_stor60: float = 0
    conv_stor61: float = 0
    conv_stor62: float = 0
    conv_stor63: float = 0
    conv_stor64: float = 0
    conv_stor65: float = 0
    conv_stor66: float = 0
    conv_stor67: float = 0
    conv_stor68: float = 0
    conv_stor69: float = 0
    conv_stor70: float = 0
    conv_stor71: float = 0
    conv_stor72: float = 0
    conv_stor73: float = 0
    conv_stor74: float = 0
    conv_stor75: float = 0
    conv_stor76: float = 0
    conv_stor77: float = 0
    conv_stor78: float = 0
    conv_stor79: float = 0
    conv_stor80: float = 0
    conv_stor81: float = 0
    conv_stor82: float = 0
    conv_stor83: float = 0
    conv_stor84: float = 0
    conv_stor85: float = 0
    conv_stor86: float = 0
    conv_stor87: float = 0
    conv_stor88: float = 0
    conv_stor89: float = 0
    conv_stor90: float = 0
    conv_stor91: float = 0
    conv_stor92: float = 0
    conv_stor93: float = 0
    conv_stor94: float = 0
    conv_stor95: float = 0
    conv_stor96: float = 0
    conv_stor97: float = 0
    conv_stor98: float = 0
    conv_stor99: float = 0


class BasinStateVariables(NamedTuple):
    index: int = 1
    name: str = "watershed"
    channelstorage: float = 0
    rivuletstorage: float = 0
    qout: tuple = (0, 0)
    qlat: tuple = (0, 0, 0, 0)
    qin: tuple = 20 * (0, )

