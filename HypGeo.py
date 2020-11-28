import sympy as sp
from typing import List, Dict, Union, Any
import matplotlib.pyplot as plt
import matplotlib.patches as mplp
import matplotlib.collections as mplc
from HypTypes import *


def generate_route_descs_2D(points: List[sp.Point2D], init_tangent: float) -> List[Dict[str, Any]]:
    """
    Generate descriptions for route arcs/segments
    Params:                                                                     type:
    :param points: Consecutive route points                                     List[sp.Point2D]
    :param init_tangent: initial tangent in first route point                   float
    :return: List of route arcs/segments descriptions                           List[Dict[str, Any]]
    """
    descs = []
    ray_tangent = sp.Ray2D(points[0], angle=init_tangent)

    for p_id in range(len(points) - 1):
        point_A = points[p_id]
        point_B = points[p_id + 1]
        arc_desc = generate_arc_2D(point_A=point_A, point_B=point_B, ray_tangent_A=ray_tangent)
        ray_tangent = arc_desc['ray_tangent_B']
        descs.append(arc_desc)

    return descs


def generate_arc_2D(point_A: sp.Point2D, point_B: sp.Point2D, ray_tangent_A: sp.Ray2D) -> Dict[str, Any]:
    """
        Function generating all parameters needed to draw an arc and tangent for next arc starting point
        Params:                                                                     type:
        :param point_A: start point                                                 Point2D
        :param point_B: end point                                                   Point2D
        :param ray_tangent_A: tangent in point A                                    Ray2D
        :return: Dictionary containing:                                             Dict[str, Any]
                    point_cR: circle center,                                        Union[Point2D, None]
                    circle_radius_len: radius length,                               Union[float, None]
                    ray_RA: ray from circle center raising to point A,              Union[Ray2D, None]
                    arc_rad_len: angle length of the arc in radians,                Union[float, None]
                    ray_tangent_B: tangent in point B,                              Ray2D]
                    segment: segment AB,                                            Union[Segment2D, None]
    """

    """ Check if AB is curve arc or segment """
    segment_AB = sp.Segment2D(point_A, point_B)
    if ray_tangent_A.is_parallel(segment_AB):
        """ AB is segment """

        """ Check if point B lies on rising direction of arc_tangent_A - if not then connection AB is impossible"""
        if not ray_tangent_A.contains(point_B):
            raise ValueError('Forbidden combination! Point B lies in straight line behind point A. '
                             'Impossible to connect points! Point A: {}, Point B: {}'.format(point_A, point_B))

        """ Evaluate tangent in point B """
        ray_tangent_B = sp.Ray2D(point_B, angle=dir_to_rad(ray_tangent_A.direction))

        return {'point_cR': sp.Point2D(np.inf, np.inf),
                'circle_radius_len': np.inf,
                'ray_RA': np.inf,
                'arc_rad_len': np.inf,
                'ray_tangent_B': ray_tangent_B,
                'segment': segment_AB}

    else:
        """ AB is arc """

        """ Evaluate circle center """
        point_cAB = segment_AB.midpoint
        line_cABR = segment_AB.perpendicular_line(point_cAB)
        line_AR = ray_tangent_A.perpendicular_line(point_A)
        point_cR = line_AR.intersection(line_cABR)[0]

        """ Evaluate length of circle radius """
        circle_radius_len = point_A.distance(point_cR)

        """ Evaluate ray starting in circle center and rising in point A direction """
        ray_RA = sp.Ray2D(point_cR, point_A)

        """ Evaluate angle length of the arc """
        ray_tangA_rot_p90 = ray_tangent_A.rotate(angle=sp.pi/2)
        ray_tangA_rot_m90 = ray_tangent_A.rotate(angle=-sp.pi/2)
        ray_AB = sp.Ray2D(point_A, point_B)

        tangA_m90_limit = dir_to_rad(direction_point=ray_tangA_rot_m90.direction)
        tangA_p90_limit = dir_to_rad(direction_point=ray_tangA_rot_p90.direction)
        dir_angle_AB = dir_to_rad(direction_point=ray_AB.direction)

        """ Check if point_cAB and point_cR are the same point - if yes then angle_ARB is pi """
        if point_cAB.equals(point_cR):
            angle_ARB = sp.pi
        else:
            angle_ARB = (sp.Triangle(point_A, point_cR, point_B)).angles[point_cR]

        if tangA_m90_limit < tangA_p90_limit:
            if tangA_m90_limit < dir_angle_AB < tangA_p90_limit:
                """ Short angle """
                arc_rad_len = angle_ARB
            else:
                """ Long angle """
                arc_rad_len = 2 * sp.pi - angle_ARB
        else:
            if dir_angle_AB < tangA_p90_limit or dir_angle_AB > tangA_m90_limit:
                """ Short angle """
                arc_rad_len = angle_ARB
            else:
                """ Long angle """
                arc_rad_len = 2 * sp.pi - angle_ARB

        """ Check rotation direction """
        if (((ray_AB.closing_angle(ray_tangent_A) < 0) and (abs(ray_AB.closing_angle(ray_tangent_A)) < sp.pi)) or
                (((ray_AB.closing_angle(ray_tangent_A)) > 0) and (abs(ray_AB.closing_angle(ray_tangent_A)) > sp.pi))):
            """ Clockwise direction - negative angle value """
            arc_rad_len = -arc_rad_len
        else:
            """ Counter clockwise direction - do nothing """
            pass

        """ Evaluate tangent in point B """
        dir_angle_tangent_A = dir_to_rad(direction_point=ray_tangent_A.direction)
        dir_angle_tangent_B = dir_angle_tangent_A + arc_rad_len
        ray_tangent_B = sp.Ray2D(point_B, angle=dir_angle_tangent_B)

        return {'point_cR': point_cR,
                'circle_radius_len': circle_radius_len,
                'ray_RA': ray_RA,
                'arc_rad_len': arc_rad_len,
                'ray_tangent_B': ray_tangent_B,
                'segment': None}


def get_route_len_2D(points: List[sp.Point2D], init_tangent: float) -> float:
    """
    Calculate route length defined by parameter points ana tangent in first point
    Params:                                                                     type:
    :param points: consecutive route points                                     List[Point2D]
    :param init_tangent: Initial tangent angle in point A in radians            float
    :return: route length                                                       float
    """
    route_len = 0
    descs = generate_route_descs_2D(points=points, init_tangent=init_tangent)

    for desc in descs:
        if desc['segment'] is not None:
            route_len += float(desc['segment'].length)
        else:
            route_len += float(2 * sp.pi * float(desc['circle_radius_len']) * (abs(float(desc['arc_rad_len'] / (2 * sp.pi)))))

    return float(route_len)


def plot_route_2d(plane: Plane, p_dicts: List[Dict[str, int]], init_tangent=sp.pi/2) -> None:
    """
    Plot route in 2D
    Params:                                                                     type:
    :param plane: Horizontal or Vertical                                        Plane
    :param p_dicts: list of consecutive route points                            List[Dict[str, int]]
    :param init_tangent: Initial tangent angle in point A in radians            float
    :return: None
    """

    ax = plt.gca()

    """ Set map limits """
    if plane == Plane.HORIZONTAL:
        ax.set_xlim(MAP_LIMIT['xmin'], MAP_LIMIT['xmax'])
        ax.set_ylim(MAP_LIMIT['ymin'], MAP_LIMIT['ymax'])
    elif plane == Plane.VERTICAL:
        route_len = (p_dicts[-1])['d'] - (p_dicts[-0])['d']
        ax.set_xlim((p_dicts[0])['d'] - 0.1 * route_len, (p_dicts[-1])['d'] + 0.1 * route_len)
        ax.set_ylim((p_dicts[0])['d'] - 0.1 * route_len, (p_dicts[-1])['d'] + 0.1 * route_len)

    """ Get route description """
    points = dict_to_points(p_dicts=p_dicts, plane=plane)
    descs = generate_route_descs_2D(points=points, init_tangent=init_tangent)

    for desc in descs:
        if desc['segment'] is None:
            """ Add arc to plot """
            xR, yR = desc['point_cR'].coordinates
            theta1 = np.degrees(float(dir_to_rad(desc['ray_RA'].direction)))
            theta2 = theta1 + np.degrees(float(desc['arc_rad_len']))
            if desc['arc_rad_len'] < 0:
                theta1, theta2 = theta2, theta1

            ax.add_patch(mplp.Arc((xR, yR), 2 * float(desc['circle_radius_len']),
                                  float(2 * desc['circle_radius_len']), theta1=float(theta1),
                                  theta2=float(theta2), edgecolor='b', lw=1.5))

            """ Add rays to plot """
            pR, pA = desc['ray_RA'].points
            pRx, pRy = pR.coordinates
            pAx, pAy = pA.coordinates
            pBx, pBy = desc['ray_tangent_B'].source.coordinates
            segment = mplc.LineCollection([[(pRx, pRy), (pAx, pAy)], [(pRx, pRy), (pBx, pBy)]], linewidths=1)
            ax.add_collection(segment)

            """ Add arc center point to plot """
            ax.scatter(pRx, pRy, cmap='viridis', linewidth=1)
            ax.text(pRx, pRy, 'r' + str(descs.index(desc)) + str(descs.index(desc) + 1), color='black')

        else:
            """ Add segment to plot """
            p1, p2 = desc['segment'].points
            p1x, p1y = p1.coordinates
            p2x, p2y = p2.coordinates
            segment = mplc.LineCollection([[(p1x, p1y), (p2x, p2y)]], linewidths=1.5)
            ax.add_collection(segment)

    """ Add points to plot """
    for point in p_dicts:
        x, y = (point['x'], point['y']) if (plane == Plane.HORIZONTAL) else (point['d'], point['z'])
        ax.scatter(x, y, cmap='viridis', linewidth=2)
        if point == p_dicts[0]:
            ax.text(x, y, 'KRK', color='black')
        elif point == p_dicts[-1]:
            ax.text(x, y, 'WAW', color='black')
        else:
            ax.text(x, y, str(p_dicts.index(point)), color='black')

    ax.grid()
    plt.show()


def dir_to_rad(direction_point: sp.Point2D) -> Union[float, sp.atan]:
    """
    Convert point indicating direction to angle in radians
    Params:                                                                     type:
    :param direction_point: indicates Ray direction in relation to point (0,0)  Point
    :return: Angle in radians                                                   Union[float, atan]
    """
    x, y = direction_point.coordinates

    if x > 0 and y > 0:
        return sp.atan(y / x)
    elif x < 0 and y > 0:
        return sp.pi + sp.atan(y / x)
    elif x < 0 and y < 0:
        return sp.pi + sp.atan(y / x)
    elif x > 0 and y < 0:
        return 2 * sp.pi + sp.atan(y / x)
    elif x > 0 and y == 0:
        return 0
    elif x == 0 and y > 0:
        return sp.pi / 2
    elif x < 0 and y == 0:
        return sp.pi
    elif x == 0 and y < 0:
        return 3 / 2 * sp.pi
    else:
        raise ValueError('Invalid point coordinates: ({},{})!'.format(x, y))


def dict_to_points(p_dicts: List[Dict[str, int]], plane: Plane) -> List[sp.Point2D]:
    """
    Convert dictionary point representation to Point2d obects
    Params:                                                                     type:
    :param p_dicts:                                                             List[Dict[str, int]]
    :param plane:                                                               Plane
    :return: List of Point2D objects                                            List[sp.Point2D]
    """
    p_out = []
    for p_dict in p_dicts:
        if plane == Plane.HORIZONTAL:
            p_out.append(sp.Point2D(p_dict['x'], p_dict['y']))
        else:
            p_out.append(sp.Point2D(p_dict['d'], p_dict['z']))
    return p_out


def evaluate_route_len(plane: Plane, p_dicts: List[Dict[str, int]], init_tangent=sp.pi/2) -> float:
    """
    Compute route length in horizontal plane
    Params:                                                                     type:
    :param plane:                                                               Plane
    :param p_dicts:                                                             List[Dict[str, int]]
    :param init_tangent:                                                        float
    :return:                                                                    float
    """
    points = dict_to_points(p_dicts=p_dicts, plane=plane)
    route_len = get_route_len_2D(points=points, init_tangent=(sp.pi/2 if plane == Plane.HORIZONTAL else 0))
    return route_len

