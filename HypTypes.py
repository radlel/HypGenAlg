from enum import Enum
import numpy as np


""" Coordinates for Krakow and Warsaw """
X_KRK = 244000
Y_KRK = 567000
Z_KRK = 208
X_WAW = 487000
Y_WAW = 638000
Z_WAW = 108

X_MIN_ALLOWED = 194000
X_MAX_ALLOWED = 536000
Y_MIN_ALLOWED = 467000
Y_MAX_ALLOWED = 809000


""" Definitins of distances """
DIST100M = 100
DIST200M = 200
DIST500M = 500
DIST1KM = 1000
DIST10KM = 10000
DIST50KM = 50000
DIST100KM = 100000
DIST1000KM = 1000000

""" Definitions of min and max route vertical level """
ROUTE_MIN_HIGHT = 0
ROUTE_MAX_HIGHT = 500

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

ANGLE_MIN_H = 1 / 12 * np.pi
ANGLE_MAX_H = 5 / 12 * np.pi
ANGLE_MIN_V = - 1 / 16 * np.pi
ANGLE_MAX_V = 1 / 16 * np.pi


""" Definitions for plotting """
PLOT_INIT = False
PLOT_FITNESS = False
PLOT_SAVE = True


PENALTY_TIGHTARC = 100
PENALTY_INVALIDROUTE = 1000000

COST_TUBE = 27 / 1000
COST_TUNNEL_BASE = 26 / 1000
COST_PYLON_PARAM = 0.094 / 1000
COST_EXC = 1 / 1000
COST_EMB = 1 / 1000


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

MAX_COST = 9999999999999999999999999999999999999999999999999999

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

