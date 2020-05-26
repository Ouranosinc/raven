"""
State variables
===============

* HRUStateVariables
* BasinStateVariables

Implemented using typing.NamedTuple

Use _replace to update individual values.

"""
from typing import NamedTuple


class HRUStateVariables(NamedTuple):
    """Initial condition for a given HRU."""
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
    """Initial conditions for a flow segment."""
    index: int = 1
    name: str = "watershed"
    channelstorage: float = 0
    rivuletstorage: float = 0
    qout: tuple = (0,)
    qoutlast: float = 0
    qlat: tuple = (0, 0, 0)
    qlatlast: float = 0
    qin: tuple = 20 * (0,)
