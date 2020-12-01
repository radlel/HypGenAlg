from enum import Enum
import numpy as np


""" Coordinates for Krakow and Warsaw """
X_KRK = 244000
Y_KRK = 567000
Z_KRK = 208
X_WAW = 487000
Y_WAW = 638000
Z_WAW = 108


""" Definitins of distances """
DIST100M = 100
DIST200M = 200
DIST500M = 500
DIST1KM = 1000
DIST10KM = 10000
DIST50KM = 50000
DIST100KM = 100000

""" Definitions of min and max route vertical level """
ROUTE_MIN_HIGHT = -500
ROUTE_MAX_HIGHT = 3000

""" Genetic Algorithm parameters """
POPULATION_SIZE = 10
CHROMOSOME_SIZE = 5
NEWPOP_BEST_PARENTS_NUM = 3
NEWPOP_BEST_PARENTS_START = 0
NEWPOP_CHILDREN_NUM = 5
NEWPOP_CHILDREN_START = 3
NEWPOP_RANDOM_NUM = 2
NEWPOP_RANDOM_START = 8
GENERATIONS_NUM = 10

""" Definition of start and end point """
START_POINT = {
    'x': X_KRK,
    'y': Y_KRK,
    'z': Z_KRK
}

END_POINT = {
    'x': X_WAW,
    'y': Y_WAW,
    'z': Z_WAW
}

""" Definitions for plotting """
PLOT_INIT = False
# PLOT_INIT = True

""" Definition of map borders """
MAP_AREA_SIZE = 343 * DIST1KM

MAP_LIMIT = {
    'xmin': X_KRK - DIST50KM,
    'xmax': X_KRK - DIST50KM + MAP_AREA_SIZE,
    'ymin': Y_KRK - DIST100KM,
    'ymax': Y_KRK - DIST100KM + MAP_AREA_SIZE,
    'zmin': -1000,
    'zmax': 2000,
    'dmin': None,
    'dmax': None
}


class Plane(Enum):
    NONE = 0,
    HORIZONTAL = 1,
    VERTICAL = 2,
    """ ========== """
    PLANE_NUM = 3

class ArcDirection(Enum):
    STRIGHT = 0,
    CLOCKWISE = 1,
    ANTICLOCKWISE = 2,
    """ ========== """
    ARCDIRECTION_NUM = 3


POINT_DEF = {
    Plane.HORIZONTAL: {'x': np.inf, 'y': np.inf},
    Plane.VERTICAL: {'z': np.inf, 'd': np.inf}
}

CROSS_SELECT = {
    0: (0, 1),
    1: (0, 2),
    2: (0, 3),
    3: (0, 4),
    4: (0, 5),
    5: (0, 6),
    6: (0, 7),
    7: (0, 8),
    8: (0, 9),
    9: (0, 10),
    10: (1, 2),
    11: (1, 3),
    12: (1, 4),
    13: (1, 5),
    14: (1, 6),
    15: (1, 7),
    16: (1, 8),
    17: (1, 9),
    18: (1, 10),
    19: (2, 3),
    20: (2, 4),
    21: (2, 5),
    22: (2, 6),
    23: (2, 7),
    24: (2, 8),
    25: (2, 9),
    26: (2, 10),
    27: (3, 4),
    28: (3, 5),
    29: (3, 6),
    30: (3, 7),
    31: (3, 8),
    32: (3, 9),
    33: (3, 10),
    34: (4, 5),
    35: (4, 6),
    36: (4, 7),
    37: (4, 8),
    38: (4, 9),
    39: (4, 10),
    40: (5, 6),
    41: (5, 7),
    42: (5, 8),
    43: (5, 9),
    44: (5, 10),
    45: (6, 7),
    46: (6, 8),
    47: (6, 9),
    48: (6, 10),
    49: (7, 8)
}

CROSS_SELECT5 = {
    0: (0, 1),
    1: (0, 2),
    2: (0, 3),
    3: (1, 2),
    4: (1, 3)
}
