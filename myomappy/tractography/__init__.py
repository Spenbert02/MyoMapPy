from enum import IntEnum


class TractographyMethods(IntEnum):
    """
    enumerations for tractography methods
    """
    EULER = 0
    RUNGE_KUTTA_4 = 1


# hard-coded maximum streamline length MAX_SL_POINTS (needs to match MAX_SL_POINTS in Streamline.h)
MAX_SL_POINTS = 1000


# import from individual tractography technique libraries
from .rungekutta import *
from .euler import *
