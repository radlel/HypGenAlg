from HypTypes import *
from copy import copy, deepcopy
from typing import Dict, Tuple, List, Union, Any
from HypGeo import get_route_description, plot_route_2d, dir_to_rad
from HypData import read_invalid_data, read_filtered_data, get_axis_z_value
import sympy as sp
import pandas as pd
import math
import matplotlib.pyplot as plt


import logging

# logging.basicConfig(filename='results.log')
logging.basicConfig(filename='output\\results.log', filemode='w', format='%(message)s')

individual_unique_id = 0
route_summary = None


""" Genetic Algorithm parameters """
CHROMOSOME_SIZE = 6
GENERATIONS_NUM = 100

NEWPOP_BEST_PARENTS_NUM = 30
NEWPOP_CHILDREN_NUM = 50
NEWPOP_RANDOM_NUM = 20

POPULATION_SIZE = NEWPOP_BEST_PARENTS_NUM + NEWPOP_CHILDREN_NUM + NEWPOP_RANDOM_NUM

NEWPOP_BEST_PARENTS_START = 0
NEWPOP_CHILDREN_START = NEWPOP_BEST_PARENTS_NUM
NEWPOP_RANDOM_START = NEWPOP_BEST_PARENTS_NUM + NEWPOP_CHILDREN_NUM

ROUTE_RESOLUTION = 1000
MATING_POINTS_MAX = 100

X_DATA_FIL, Y_DATA_FIL, Z_DATA_FIL = np.inf, np.inf, np.inf
X_DATA_INV, Y_DATA_INV = np.inf, np.inf


class Gen:
    """ Definition of Gen - coordinates of single point """
    def __init__(self, plane: Plane) -> None:
        self.point = deepcopy(POINT_DEF[plane])


class Chromosome:
    """ Definition of Chromosome - collection of Gens """
    def __init__(self, plane: Plane) -> None:
        self.plane = plane
        self.gens = np.array([Gen(plane=plane) for i in range(CHROMOSOME_SIZE)])
        self.gen_init_tangent = np.inf

    def init_tangent_rand(self):
        if self.plane == Plane.HORIZONTAL:
            """ Initialize tangent in first point (for horizontal drawn from range (ANGLE_MIN_H, ANGLE_MAX_H) """
            self.gen_init_tangent = (ANGLE_MAX_H - ANGLE_MIN_H) * np.random.sample() + ANGLE_MIN_H
        elif self.plane == Plane.VERTICAL:
            """ Initialize tangent in first point (for horizontal drawn from range (ANGLE_MIN_V, ANGLE_MAX_V)) """
            self.gen_init_tangent = (ANGLE_MAX_V - ANGLE_MIN_V) * np.random.sample() + ANGLE_MIN_V
        else:
            raise ValueError('Invalid plane parameter value! {}'.format(self.plane))


class Genotype:
    """ Definition of Genotype - collection of Chromosomes """
    def __init__(self) -> None:
        self.chromosome_h = Chromosome(plane=Plane.HORIZONTAL)
        self.checksum_h = np.inf
        self.chromosome_v = Chromosome(plane=Plane.VERTICAL)
        self.checksum_v = np.inf

    def init_horizontal_random(self) -> None:
        """
        Initialize Gens (Points) for hotizontal movement
        Params:                                                                     type:
        :return: None                                                               None
        """

        """ Initialize start point and end point coordinates """
        gens = self.chromosome_h.gens
        gens[0].point['x'] = START_POINT['x']
        gens[0].point['y'] = START_POINT['y']
        gens[-1].point['x'] = END_POINT['x']
        gens[-1].point['y'] = END_POINT['y']

        """ Draw intermediate points coordinates """
        for gen_id in range(1, CHROMOSOME_SIZE - 1):
            while True:
                # x_drawn = np.random.choice([i for i in range(MAP_LIMIT['xmin'] + DIST10KM,
                #                                              MAP_LIMIT['xmax'] - DIST10KM)], 1)[0]
                # y_drawn = np.random.choice([i for i in range(MAP_LIMIT['ymin'] + DIST10KM,
                #                                              MAP_LIMIT['ymax'] - DIST10KM)], 1)[0]

                x_drawn = np.random.choice([i for i in range(X_KRK, X_WAW)], 1)[0]
                y_drawn = np.random.choice([i for i in range(Y_KRK, Y_WAW)], 1)[0]

                """ Check if selected point has proper z value (not np.inf) """
                if is_point_valid(x=x_drawn, y=y_drawn):
                    gens[gen_id].point['x'] = x_drawn
                    gens[gen_id].point['y'] = y_drawn
                    break
                else:
                    """" Drawn point not in range - draw another one """
                    pass

        """ For better convergence at the beginning sort points to be from most far away to the closest end point """
        values_x = deepcopy([gen.point['x'] for gen in gens[1:-1]])
        values_y = deepcopy([gen.point['y'] for gen in gens[1:-1]])
        values_x.sort()
        values_y.sort()

        for gen_id in range(1, len(gens) - 1):
            (gens[gen_id]).point['x'] = values_x[gen_id - 1]
            (gens[gen_id]).point['y'] = values_y[gen_id - 1]

        """ Initialize tangent in first point (for horizontal drawn from range (0, pi/2))"""
        self.chromosome_h.init_tangent_rand()

        """ Calculate checksum based on points and tangent """
        self.init_checksum(plane=Plane.HORIZONTAL)

    def init_vertical_random(self) -> None:
        """
        Initialize Gens (Points) for vertical movement
        Params:                                                                     type:
        :return: None
        """

        """ Initialize start point and end point coordinates """
        gens = self.chromosome_v.gens
        gens[0].point['z'] = START_POINT['z']
        gens[0].point['d'] = 0
        gens[-1].point['z'] = END_POINT['z']
        gens[-1].point['d'] = ROUTE_RESOLUTION

        """ Draw intermediate points coordinates; z - height in adequate point d,
                                                  d - point in route len divided for 1000 equal segemnts """
        z_drawns = np.random.choice([i for i in range(ROUTE_MIN_HIGHT, ROUTE_MAX_HIGHT)], CHROMOSOME_SIZE - 2)

        """ Make sure d values are unique """
        while True:
            d_drawns = np.random.choice(ROUTE_RESOLUTION, CHROMOSOME_SIZE - 2)

            if len(d_drawns) == len(set(d_drawns)):
                """" d values are unique - ok """
                break
            else:
                """ d values are not unique, draw another collection """
                pass

        """ Sort values by growing distance d """
        df_vert = pd.DataFrame({'z': z_drawns, 'd': d_drawns})
        df_vert = deepcopy(df_vert.sort_values(by=['d'], ignore_index=True))

        """ Set intermediate points """
        for gen_id in range(1, CHROMOSOME_SIZE - 1):
            gens[gen_id].point['z'] = df_vert['z'][gen_id - 1]
            gens[gen_id].point['d'] = df_vert['d'][gen_id - 1]

        """ Initialize tangent in first point """
        self.chromosome_v.init_tangent_rand()

        """ Calculate checksum based on points and tangent """
        self.init_checksum(plane=Plane.VERTICAL)

    def get_points_descs(self, plane: Plane) -> List[Dict[str, int]]:
        if plane == Plane.HORIZONTAL:
            return [{'x': gen.point['x'], 'y': gen.point['y']} for gen in self.chromosome_h.gens]
        elif plane == Plane.VERTICAL:
            return [{'z': gen.point['z'], 'd': gen.point['d']} for gen in self.chromosome_v.gens]
        else:
            raise ValueError('Invalid plane parameter value! {}'.format(plane))

    def init_checksum(self, plane: Plane):
        """ Calculate checksum - only intermediate points plus tangent value times 100 """
        if plane == Plane.HORIZONTAL:
            self.checksum_h = (math.floor(sum([float(gen.point['x'] + gen.point['y']) for gen in
                               self.chromosome_h.gens[1:-1]]) + float(self.chromosome_h.gen_init_tangent) * 100))
            return copy(self.checksum_h)
        elif plane == Plane.VERTICAL:
            self.checksum_v = (math.floor(sum([float(gen.point['z'] + gen.point['d']) for gen in
                               self.chromosome_v.gens[1:-1]]) + float(self.chromosome_v.gen_init_tangent) * 100))
            return copy(self.checksum_v)


class ArcDesc:
    def __init__(self):
        self.point_cR = sp.Point2D(np.inf, np.inf)
        self.circle_radius_len = np.inf
        self.ray_RA = sp.Ray2D(sp.Point2D(0, 0), angle=0)
        self.arc_rad_len = np.inf
        self.ray_tangent_B = sp.Ray2D(sp.Point2D(0, 0), angle=0)
        self.segment = np.inf


class RouteDesc:
    def __init__(self):
        self.arcs = np.array([ArcDesc() for i in range(CHROMOSOME_SIZE - 1)])


class Fenotype:
    def __init__(self):
        self.route_horizontal = RouteDesc()
        self.route_vertical = RouteDesc()
        self.route_len_h = np.inf
        self.route_len_total = np.inf
        self.route_desc_h = np.inf
        self.route_desc_v = np.inf
        self.fitness_val = np.inf

    def init_horizontal(self, p_dicts: List[Dict[str, int]], init_tangent: float) -> None:
        route_desc, route_len = get_route_description(plane=Plane.HORIZONTAL, p_dicts=p_dicts, init_tangent=init_tangent)
        arcs = self.route_horizontal.arcs
        for (arc, arc_desc) in zip(arcs, route_desc):
            arc.point_cR = arc_desc['point_cR']
            arc.circle_radius_len = arc_desc['circle_radius_len']
            arc.ray_RA = arc_desc['ray_RA']
            arc.arc_rad_len = arc_desc['arc_rad_len']
            arc.ray_tangent_B = arc_desc['ray_tangent_B']
            arc.segment = arc_desc['segment']

        self.route_len_h = route_len
        self.route_desc_h = route_desc

        if PLOT_INIT:
            global individual_unique_id
            plot_route_2d(plane=Plane.HORIZONTAL, route_desc=route_desc, route_len=route_len, p_dicts=p_dicts,
                          title=str(individual_unique_id) + 'h')

    def init_vertical(self, p_dicts: List[Dict[str, int]], init_tangent: float) -> None:
        route_desc, route_len = get_route_description(plane=Plane.VERTICAL, p_dicts=p_dicts, init_tangent=init_tangent)
        arcs = self.route_vertical.arcs
        for (arc, arc_desc) in zip(arcs, route_desc):
            arc.point_cR = arc_desc['point_cR']
            arc.ray_RA = arc_desc['ray_RA']
            arc.arc_rad_len = arc_desc['arc_rad_len']
            arc.ray_tangent_B = arc_desc['ray_tangent_B']
            arc.segment = arc_desc['segment']

        self.route_len_total = route_len
        self.route_desc_v = route_desc

        if PLOT_INIT:
            global individual_unique_id
            plot_route_2d(plane=Plane.VERTICAL, route_desc=route_desc, route_len=route_len, p_dicts=p_dicts,
                          title=str(individual_unique_id) + 'v')
            individual_unique_id += 1

    def get_route_desc(self, plane: Plane):
        if plane == Plane.HORIZONTAL:
            return self.route_desc_h
        elif plane == Plane.VERTICAL:
            return self.route_desc_v
        else:
            raise ValueError('Invalid plane parameter value! {}'.format(plane))


class Individual:
    def __init__(self):
        self.genotype = Genotype()
        self.fenotype = Fenotype()
        self.discete_route_h_points = None

    def init_fenotype_is_valid(self) -> bool:
        """ Assumed genotypes both h and v are already initialized at this point
            In here only initialize fenotypes and check if Individual is valid """

        self.fenotype.init_horizontal(p_dicts=self.genotype.get_points_descs(plane=Plane.HORIZONTAL),
                                      init_tangent=self.genotype.chromosome_h.gen_init_tangent)
        disc_route_ret = self.dicretize_route()

        """ Check if all discrete route points are in known region """
        if disc_route_ret is not None:
            """ Route points are in known region - ok """
            self.discete_route_h_points = disc_route_ret
        else:
            """ Route points not in region - return False """
            print('- ' * 8 + 'Individual: Route in unknown region, Individual invalid!')
            return False

        self.fenotype.init_vertical(p_dicts=self.map_genotype_v_to_p_dicts(),
                                    init_tangent=self.genotype.chromosome_v.gen_init_tangent)

        for arc_desc in self.fenotype.route_desc_v:
            if arc_desc['arc_rad_len'] > np.pi:
                print('- ' * 8 + 'Individual: arc angle too big: {}, Individual invalid!'
                      .format(np.degrees(float(arc_desc['arc_rad_len']))))
                """ Free discrete points as they are not needed any longer and only occupy RAM """
                self.discete_route_h_points = None
                return False

        """ If fenotype is ok then calculate fitness """
        print('. ' * 8 + 'Individual: Fenotype ok, calculate fitness')
        self.calculate_fitness()

        """ Free discrete points as they are not needed any longer and only occupy RAM """
        self.discete_route_h_points = None
        return True

    def initialize_random(self) -> None:
        while True:
            self.genotype.init_horizontal_random()
            self.genotype.init_vertical_random()
            individual_valid = self.init_fenotype_is_valid()

            if individual_valid is True:
                print('. ' * 8 + 'Individual: Initialization was successfull')
                break
            else:
                print('- ' * 8 + 'Individual: Initialization was not successfull, retry')

    def map_genotype_v_to_p_dicts(self) -> List[Dict[str, int]]:
        points_v = deepcopy(self.genotype.get_points_descs(plane=Plane.VERTICAL))
        route_len_h = copy(self.fenotype.route_len_h)
        return [{'z': point['z'], 'd': math.floor(point['d'] / ROUTE_RESOLUTION * route_len_h)} for point in points_v]

    def calculate_fitness(self) -> None:

        """ Take dicretized route """
        route_points = self.discete_route_h_points

        """ Get vector of hight points from generatd route """
        gen_z_vals = self.get_gen_route_hights(route_points=route_points, route_v_descs=self.fenotype.route_desc_v)

        """ Get vector of hight points from landform """
        orig_landform_z_vals = self.get_landform_route_heights(route_points=route_points)

        if len(gen_z_vals) != len(orig_landform_z_vals):
            raise ValueError('Len of arrays must be equal!')

        """ Create vector of hight differences """
        diff_vector_z = [gen_z_vals[i] - orig_landform_z_vals[i] for i in range(len(gen_z_vals))]

        cost, summary = self.calculate_route_cost(diff_vector_z=diff_vector_z)

        """ Check if arc is not to tight, in case it is add penalty """
        radius_lens = [int(arc['circle_radius_len']) for arc in self.fenotype.route_desc_h]
        for radius_len in radius_lens:
            if radius_len < 23000:
                print('DEBUG Fitness: Arc to tight, add PENALTY_TIGHTARC, radius:', radius_lens)
                cost += PENALTY_TIGHTARC
                break

        print('. ' * 8 + 'Individual: Calculated route cost:', format_fitness(fitness=cost))

        self.fenotype.fitness_val = cost

        if PLOT_SAVE:
            global individual_unique_id
            plot_route_2d(plane=Plane.HORIZONTAL,
                          route_desc=self.fenotype.route_desc_h,
                          route_len=self.fenotype.route_len_h,
                          p_dicts=self.genotype.get_points_descs(plane=Plane.HORIZONTAL),
                          title=str(individual_unique_id) + 'h')

            p_dicts = self.map_genotype_v_to_p_dicts()
            route_desc, route_len = get_route_description(plane=Plane.VERTICAL,
                                                          p_dicts=p_dicts,
                                                          init_tangent=self.genotype.chromosome_v.gen_init_tangent)

            d_points_disc = [elem['d'] for elem in route_points]
            landform = (d_points_disc, orig_landform_z_vals)
            z_min = min(gen_z_vals)
            z_max = max(gen_z_vals)
            plot_route_2d(plane=Plane.VERTICAL,
                          route_desc=route_desc,
                          route_len=route_len,
                          p_dicts=p_dicts,
                          title=str(individual_unique_id) + 'v',
                          landform=landform,
                          z_min=z_min,
                          z_max=z_max)


            logging.warning(str(individual_unique_id) + ' ' +
                            str(int(cost)) + ' ' +
                            str(summary) + ' ' +
                            str(self.fenotype.route_desc_h) + ' ' +
                            str(self.fenotype.route_desc_v))
            individual_unique_id += 1




    def calculate_route_cost(self, diff_vector_z: List[float]) -> (int, Dict):
        cost = 0
        tu, ex, gr, em, py = 0, 0, 0, 0, 0,

        for diff in diff_vector_z:

            """ Check if point is valid and has no np.inf valid iside """
            if diff == -np.inf:
                print('DEBUG Fitness: Invalid route, return PENALTY_INVALIDROUTE, points:', diff_vector_z)
                return MAX_COST, {'tu': 0, 'ex': 0, 'gr': 0, 'em': 0, 'py': 0}

            try:

                if diff <= -6:
                    """ tunnel """
                    cost += COST_TUBE + COST_TUNNEL_BASE * (2 ** math.floor(abs(diff) / 50))
                    tu += 1
                elif -6 < diff <= -1:
                    """ excavation """
                    cost += COST_TUBE + COST_EXC
                    ex += 1
                elif -1 < diff <= 1:
                    """ on ground """
                    cost += COST_TUBE
                    gr += 1
                elif 1 < diff <= 6:
                    """ embankment """
                    cost += COST_TUBE + COST_EMB
                    em += 1
                elif diff > 6:
                    """ pylon """
                    cost += COST_TUBE + COST_PYLON_PARAM * (diff ** 2)
                    py += 1
                else:
                    print(diff)
                    raise ValueError

                """ Check if intiger is not about to limit """
                if cost >= MAX_COST:
                    print('DEBUG Fitness: Reached MAX COST!')
                    tu, ex, gr, em, py = 0, 0, 0, 0, 0,
                    return MAX_COST, {'tu': 0, 'ex': 0, 'gr': 0, 'em': 0, 'py': 0}

            except OverflowError:
                print('DEBUG OverflowError: int too large to convert to float, return MAX COST')
                return MAX_COST, {'tu': 0, 'ex': 0, 'gr': 0, 'em': 0, 'py': 0}

        """ Add mainenance costs in 10 years """
        cost *= 2

        return cost, {'tu': tu, 'ex': ex, 'gr': gr, 'em': em, 'py': py}


    def is_route_in_region(self, route_points: List[Dict[str, Union[int, float]]]) -> bool:
        x_collection = [elem['x'] for elem in route_points]
        y_collection = [elem['y'] for elem in route_points]
        x_min = min(x_collection)
        x_max = max(x_collection)
        y_min = min(y_collection)
        y_max = max(y_collection)

        if (x_min >= X_MIN_ALLOWED and x_max <= X_MAX_ALLOWED and y_min >= Y_MIN_ALLOWED and y_max <= Y_MAX_ALLOWED and
            is_point_valid(x=x_min, y=y_min) and
            is_point_valid(x=x_min, y=y_max) and
            is_point_valid(x=x_max, y=y_max) and
            is_point_valid(x=x_max, y=y_min)):
            return True
        else:
            return False

    def get_landform_route_heights(self, route_points: List[Dict[str, Union[int, float]]]) -> List[float]:
        z_orig_vals = []

        for point in route_points:
            x, y = round(point['x'], -3), round(point['y'], -3)

            if x == X_MIN_ALLOWED - DIST1KM:
                x = X_MIN_ALLOWED
            elif x == X_MAX_ALLOWED + DIST1KM:
                x = X_MAX_ALLOWED

            if y == Y_MIN_ALLOWED - DIST1KM:
                y = Y_MIN_ALLOWED
            elif y == Y_MAX_ALLOWED + DIST1KM:
                y = Y_MAX_ALLOWED

            z_orig_vals.append(get_axis_z_value(z_data=Z_DATA_FIL,
                                                x_coordinate=x,
                                                y_coordinate=y,
                                                x_req=X_DATA_FIL,
                                                y_req=Y_DATA_FIL))

        return z_orig_vals

    def get_gen_route_hights(self, route_points: List[Dict[str, Union[int, float]]], route_v_descs: List[Dict[str, Any]]) ->\
            List[float]:

        """ Evaluate common points of fallowing arcs """
        d_points_mapped = [d_elem['d'] for d_elem in self.map_genotype_v_to_p_dicts()]

        generated_route_z_values = []

        for d_point in [math.floor(d_elem['d']) for d_elem in route_points]:
            if d_points_mapped[0] <= d_point < d_points_mapped[1]:
                """ Take first arc for calcultaions """
                arc_id = 0
                pass
            elif d_points_mapped[1] <= d_point < d_points_mapped[2]:
                """ Take second arc for calcultaions """
                arc_id = 1
            elif d_points_mapped[2] <= d_point < d_points_mapped[3]:
                """ Take third arc for calcultaions """
                arc_id = 2
            elif d_points_mapped[3] <= d_point <= d_points_mapped[4]:
                """ Take fourth arc for calcultaions """
                arc_id = 3
            elif d_points_mapped[4] <= d_point <= d_points_mapped[5]:
                arc_id = 4
            elif d_points_mapped[5] <= d_point <= d_points_mapped[6]:
                arc_id = 5
            else:
                print(d_points_mapped)
                print(d_point)
                raise ValueError('Only 5 points considered in here!')

            z = self.get_gen_disc_point_hight(arc_desc=route_v_descs[arc_id], d_point=d_point)
            generated_route_z_values.append({'d': d_point, 'z': z})

        if PLOT_FITNESS:
            plot_route_2d(plane=Plane.VERTICAL,
                          route_desc=route_v_descs,
                          route_len=self.fenotype.route_len_total,
                          p_dicts=generated_route_z_values)

        return [elem['z'] for elem in generated_route_z_values]

    def get_gen_disc_point_hight(self, arc_desc: Dict[str, Any], d_point: float):
        d_center, z_center = arc_desc['point_cR'].coordinates
        direction = ArcDirection.ANTICLOCKWISE if arc_desc['arc_rad_len'] >= 0 else ArcDirection.CLOCKWISE

        alpha = np.arccos(float((d_point - d_center) / arc_desc['circle_radius_len']))
        if direction == ArcDirection.ANTICLOCKWISE:
            alpha += np.pi

        return np.int(arc_desc['circle_radius_len'] * np.sin(alpha) + z_center)

    def dicretize_route(self) -> Union[List[Dict[str, Union[int, float]]], None]:
        curr_distance = 0.0
        dicrete_route_points = []

        for arc_desc in self.fenotype.route_desc_h:
            curr_distance, points = self.discretize_arc(arc_desc=arc_desc, curr_dist=curr_distance)
            dicrete_route_points += points

        dicrete_route_points.append({'x': X_WAW, 'y': Y_WAW, 'd': math.floor(self.fenotype.route_len_h)})

        in_region = self.is_route_in_region(route_points=dicrete_route_points)

        if in_region is True:
            return dicrete_route_points
        else:
            return None

    def discretize_arc(self, arc_desc: Dict[str, Any], curr_dist: float) -> Tuple[float, List]:
        arc_rad_len = float(arc_desc['arc_rad_len'])
        circle_radius_len = float(arc_desc['circle_radius_len'])
        direction = ArcDirection.ANTICLOCKWISE if arc_rad_len >= 0 else ArcDirection.CLOCKWISE
        x_center, y_center = arc_desc['point_cR'].coordinates

        points = []
        distance = curr_dist
        arc_dist = 0.0
        ang_1km_fraction = DIST1KM / float(2 * sp.pi * circle_radius_len)
        ang_1km_rad = float(2 * np.pi * ang_1km_fraction)
        base_angle = float(dir_to_rad(arc_desc['ray_RA'].direction))

        n_1km_steps = math.ceil(abs(arc_rad_len) / ang_1km_rad)

        for p_id in range(n_1km_steps):
            if direction == ArcDirection.ANTICLOCKWISE:
                p_angle = base_angle + p_id * ang_1km_rad
            else:
                p_angle = base_angle - p_id * ang_1km_rad

            y = circle_radius_len * np.sin(p_angle) + y_center
            x = circle_radius_len * np.cos(p_angle) + x_center

            points.append({'x': x, 'y': y, 'd': distance})

            if p_id != n_1km_steps - 1:
                dist_temp = float(ang_1km_fraction * 2 * np.pi * circle_radius_len)
                distance += dist_temp
                arc_dist += dist_temp
            else:
                ang_rest = abs(arc_rad_len) - (n_1km_steps - 1) * ang_1km_rad
                rest_dist = float((ang_rest / (2 * np.pi)) * 2 * np.pi * circle_radius_len)
                distance += rest_dist
                arc_dist += rest_dist

        return distance, points


class Population:
    """ Definition of Population - collection of Genotypes """
    def __init__(self, pop_size=POPULATION_SIZE) -> None:
        self.individuals = np.array([Individual() for i in range(pop_size)])

    def initialize_random(self, curr_pop_crcs: Union[Tuple[List, List], None] = None) -> None:
        """
        Initialize Individuals with random values
        Params:                                                                     type:
        :return: None                                                               None
        """
        print('Population: started random initialization')
        """ Get invalid coordinates so Genotype can check if drawn point is correct """

        for (individual, ind_id) in zip(self.individuals, range((len(self.individuals)))):
            print('Population: initialize Indywidual {}'.format((list(self.individuals)).index(individual)))
            while True:
                individual.initialize_random()

                """ Check if individual CRC is unique in context of this Population """
                if individual.genotype.checksum_h not in \
                        [ind.genotype.checksum_h for ind in self.individuals[:ind_id]] and \
                        individual.genotype.checksum_v not in \
                        [ind.genotype.checksum_v for ind in self.individuals[:ind_id]]:

                    """ If yes then check if is unique in context of super population - if exist """
                    if curr_pop_crcs is not None:
                        spr_pop_crcs_h, spr_pop_crcs_v = curr_pop_crcs
                        if individual.genotype.checksum_h not in spr_pop_crcs_h and \
                                individual.genotype.checksum_v not in spr_pop_crcs_v:
                            """ Inidividual is globally unique - init next individual """
                            break
                        else:
                            print('***Population: At least one of CRCs exists in SUPER population, retry random init')
                    else:
                        """ No super population - individual is unique """
                        break
                else:
                    print('***Population: At least one of CRCs exists in THIS population, retry random init')

        print('Population: ended random initialization')

        print_population_info('Random initialization', pop=self.individuals)

    def create_new_generation(self) -> None:
        print('\n', 'STARTNG CREATION OF NEW GENERATION...')

        """ Create new offspring """
        offspring = np.array([Individual() for i in range(len(self.individuals))])

        """ select 30% best parents """
        best_to_worst_parents = self.__get_best_individuals_sorted()

        print_population_info(title='Best parents', pop=best_to_worst_parents[:NEWPOP_BEST_PARENTS_NUM])

        offspring[NEWPOP_BEST_PARENTS_START:NEWPOP_BEST_PARENTS_NUM] =\
            deepcopy(best_to_worst_parents[0:NEWPOP_BEST_PARENTS_NUM])

        """ use best parents to cross over offspring and generate 50% of them """
        children = self.__cross_over_mutate(parents_sorted=best_to_worst_parents)

        print_population_info(title='Children after crossing over and mutation', pop=children)

        offspring[NEWPOP_CHILDREN_START:NEWPOP_CHILDREN_START + NEWPOP_CHILDREN_NUM] = deepcopy(children)

        """ create random 20% individuals """
        random_pop = Population(pop_size=NEWPOP_RANDOM_NUM)

        """ Get current offspring CRCs to make sure created randomly are not the same """
        offspr_crcs_h = [ind.genotype.checksum_h for ind in offspring]
        offspr_crcs_v = [ind.genotype.checksum_v for ind in offspring]

        random_pop.initialize_random(curr_pop_crcs=(offspr_crcs_h, offspr_crcs_v))

        offspring[NEWPOP_RANDOM_START:NEWPOP_RANDOM_START + NEWPOP_RANDOM_NUM] = deepcopy(random_pop.individuals)

        """ Set new population """
        for (individual, index) in zip(self.individuals, range(len(self.individuals))):
            individual.genotype = deepcopy(offspring[index].genotype)
            individual.fenotype = deepcopy(offspring[index].fenotype)

        print_population_info(title='New population', pop=self.individuals)

        # debug
        pop_crcs_h = [ind.genotype.checksum_h for ind in self.individuals]
        pop_crcs_v = [ind.genotype.checksum_v for ind in self.individuals]
        if len(pop_crcs_h) == len(set(pop_crcs_h)):
            """ ok """
            pass
        else:
            raise RuntimeError('H CRCs not unique!: {}'.format(pop_crcs_h))
        if len(pop_crcs_v) == len(set(pop_crcs_v)):
            """ ok """
            pass
        else:
            raise RuntimeError('V CRCs not unique!: {}'.format(pop_crcs_v))
        # end debug

        print('...FINISHED CREATION OF NEW GENERATION')

    def __get_best_individuals_sorted(self) -> np.array:
        individuals = deepcopy(self.individuals)
        costs_desc = deepcopy([{'index': index, 'cost': individual.fenotype.fitness_val} for
                              (index, individual) in zip(range(len(individuals)), individuals)])

        sorted_indexes = [index for (cost, index) in sorted(zip([item['cost'] for item in costs_desc],
                                                                [item['index'] for item in costs_desc]))]

        return np.array([individuals[index] for index in sorted_indexes])

    def get_best_individual(self) -> Individual:
        return (self.__get_best_individuals_sorted())[0]

    def __cross_over_mutate(self, parents_sorted: np.array) -> np.array:
        """ """
        print('\t\t\tStarted cross over...')
        """ Create new empty population - children """
        pop = Population(pop_size=NEWPOP_CHILDREN_NUM)
        children = pop.individuals
        parents = deepcopy(parents_sorted)
        child_id = 0

        """ Create list of mating ranges for every parent - every parent index from list parents indicates
            top value of selection range, bottom selection range is previous value in list """
        parents_mating_points = []
        for parent_id in range(len(parents)):
            parents_mating_points.append(math.ceil(parents[0].fenotype.fitness_val /
                                                   parents[parent_id].fenotype.fitness_val *
                                                   MATING_POINTS_MAX)
                                         + (parents_mating_points[parent_id - 1] if parent_id > 0 else 0))
        mating_points_sum = parents_mating_points[-1]
        print('\t\t\t\tParents mating points:', parents_mating_points, 'sum:', mating_points_sum)

        """ Process crossing procedure till all children are produced """
        while True:
            """ Draw parents by 2 random points from range (0, mating_points_sum) """
            rand_mating_points = np.random.choice(mating_points_sum, 2)
            mating_point1 = rand_mating_points[0]
            mating_point2 = rand_mating_points[1]
            parent1, parent2 = None, None
            child = Individual()

            """ Evaluate parents """
            for mate_range, parent in zip(parents_mating_points, parents):
                if mating_point1 < mate_range and parent1 is None:
                    parent1 = parent

                if mating_point2 < mate_range and parent2 is None:
                    parent2 = parent

                """ Check if parents already chosen """
                if parent1 is not None and parent2 is not None:
                    break

            print('\t\t\tDrawn parents: [{} {}]:, {} in ({},{}), {} in ({},{})'.
                  format(list(parents).index(parent1), list(parents).index(parent2), mating_point1,
                  parents_mating_points[list(parents).index(parent1) - 1],
                  parents_mating_points[list(parents).index(parent1)], mating_point2,
                  parents_mating_points[list(parents).index(parent2) - 1],
                  parents_mating_points[list(parents).index(parent2)]))

            """ Check if chosen parents are not the same Individual """
            if ((parent1.genotype.checksum_h != parent2.genotype.checksum_h) and
                    (parent1.genotype.checksum_v != parent2.genotype.checksum_v)):

                """ Generate crossing masks and make sure it is not uniform """
                while True:
                    cross_mask_h = np.random.choice([0, 1], CHROMOSOME_SIZE - 2)
                    if not sum(cross_mask_h) in [0, len(cross_mask_h)]:
                        break

                while True:
                    cross_mask_v = np.random.choice([0, 1], CHROMOSOME_SIZE - 2)
                    if not sum(cross_mask_v) in [0, len(cross_mask_v)]:
                        break

                """ Generate mask for choosing init tangent in starting point"""
                tang_mask_h = np.random.choice([0, 1], 1)[0]
                tang_mask_v = np.random.choice([0, 1], 1)[0]

                """ Start and end point is always the same """
                child.genotype.chromosome_h.gens[0] = deepcopy(parent1.genotype.chromosome_h.gens[0])
                child.genotype.chromosome_h.gens[-1] = deepcopy(parent1.genotype.chromosome_h.gens[-1])
                child.genotype.chromosome_v.gens[0] = deepcopy(parent1.genotype.chromosome_v.gens[0])
                child.genotype.chromosome_v.gens[-1] = deepcopy(parent1.genotype.chromosome_v.gens[-1])

                """ Copy intermediate gens according to crossing mask """
                for i in range(1, CHROMOSOME_SIZE - 1):
                    child.genotype.chromosome_h.gens[i] = \
                        deepcopy((parent1.genotype.chromosome_h.gens[i]) if (cross_mask_h[i - 1] == 0)
                                 else (parent2.genotype.chromosome_h.gens[i]))

                    child.genotype.chromosome_v.gens[i] = \
                        deepcopy((parent1.genotype.chromosome_v.gens[i]) if (cross_mask_v[i - 1] == 0)
                                 else (parent2.genotype.chromosome_v.gens[i]))

                """ Vertical points must be ascending - sort them if they are not """
                z_vals = [gen.point['z'] for gen in child.genotype.chromosome_v.gens[1:-1]]
                d_vals = [gen.point['d'] for gen in child.genotype.chromosome_v.gens[1:-1]]

                """ d values must be unique """
                if len(d_vals) == len(set(d_vals)):
                    """" d values are unique - ok """

                    df_vert = pd.DataFrame({'z': z_vals, 'd': d_vals})
                    df_vert = deepcopy(df_vert.sort_values(by=['d'], ignore_index=True))
                    """ Set intermediate points """
                    for gen_id in range(1, CHROMOSOME_SIZE - 1):
                        child.genotype.chromosome_v.gens[gen_id].point['z'] = deepcopy(df_vert['z'][gen_id - 1])
                        child.genotype.chromosome_v.gens[gen_id].point['d'] = deepcopy(df_vert['d'][gen_id - 1])

                    child.genotype.chromosome_h.gen_init_tangent = deepcopy(
                        parent1.genotype.chromosome_h.gen_init_tangent if (tang_mask_h == 0)
                        else parent2.genotype.chromosome_h.gen_init_tangent)

                    child.genotype.chromosome_v.gen_init_tangent = deepcopy(
                        parent1.genotype.chromosome_v.gen_init_tangent if (tang_mask_v == 0)
                        else parent2.genotype.chromosome_v.gen_init_tangent)

                    """ Process mutation """
                    child = self.__mutate(child=child)

                    """ Check if created individual does not already exist """
                    crc_h_child = child.genotype.init_checksum(plane=Plane.HORIZONTAL)
                    crc_v_child = child.genotype.init_checksum(plane=Plane.VERTICAL)

                    if ((crc_h_child not in [ind.genotype.checksum_h for ind in parents] and
                         (crc_v_child not in [ind.genotype.checksum_v for ind in parents])) and
                            (crc_h_child not in [ind.genotype.checksum_h for ind in children[:child_id]]) and
                            (crc_v_child not in [ind.genotype.checksum_v for ind in children[:child_id]])):

                        """ Initialize individual fenotype and check if is valid """
                        individual_valid = child.init_fenotype_is_valid()
                        if individual_valid is True:
                            print('. ' * 8 + 'Population: Individual valid, crossing successfull')

                            """ Add individual to new population """
                            children[child_id] = deepcopy(child)
                            print('\t\t\t\tChild crossing successfull:', child_id)
                            child_id += 1

                            if child_id == NEWPOP_CHILDREN_NUM:
                                """ Produced all children, return """
                                print('\t\t\t...Ended cross over')
                                return children

                            else:
                                """ Not all children produced yet - go to beginning of the loop """
                                pass
                        else:
                            """ Individual invalid - go to beginning of the loop """
                            print('- ' * 8 + '*** Individual has invalid fenotype, retry crossing procedure')

                    else:
                        """ Individual already exists in population - go to beginning of the loop """
                        print('\t\t\t*** Individual with at least one the same CRC already exists in population, draw '
                              'another parents pair')
                else:
                    """ Invalid d values collection, try again """
                    print('\t\t\t*** Invalid d values collection, try again')
            else:
                """ Invalid parents drawing - go to beginning of the loop """
                print('\t\t\t*** Parents has at least one same CRC, draw another parents pair')

    def __mutate(self, child: np.array) -> np.array:

        """ Horizontal mutation """
        mutate_mask = np.random.choice(10, 1)[0]

        """ Drawning 0 is 10% chance - mutations probability """
        if mutate_mask == 0:
            """ Draw which gen will be affected """

            mut_gen_id = (np.random.choice(CHROMOSOME_SIZE - 2, 1))[0] + 1

            """ Draw new point and replace """
            x_drawn = np.random.choice([i for i in range(X_KRK, X_WAW)], 1)[0]
            y_drawn = np.random.choice([i for i in range(Y_KRK, Y_WAW)], 1)[0]

            print('\t\t\tMutation will be applied for HORIZONTAL: {} by {} on place {}'
                  .format([(gen.point['x'], gen.point['y'])for gen in child.genotype.chromosome_h.gens],
                          (x_drawn, y_drawn), mut_gen_id))

            child.genotype.chromosome_h.gens[mut_gen_id].point['x'] = x_drawn
            child.genotype.chromosome_h.gens[mut_gen_id].point['y'] = y_drawn

        else:
            """ Children is not about to be mutated in horizontal, just pass """
            pass

        """ Vertical mutation """
        mutate_mask = np.random.choice(10, 1)[0]
        """ Drawning 0 is 10% chance - mutations probability """
        if mutate_mask in [0, 1, 2, 3, 4]:
            """ Draw which gen will be affected """

            mut_gen_id = (np.random.choice(CHROMOSOME_SIZE - 2, 1))[0] + 1

            """ Draw new z value and replace """
            z_drawn = np.random.choice([i for i in range(ROUTE_MIN_HIGHT, ROUTE_MAX_HIGHT)], 1)[0]

            print('\t\t\tMutation will be applied for VERTICAL: {} by {} on place {}'
                  .format([gen.point['z'] for gen in child.genotype.chromosome_v.gens],
                          z_drawn, mut_gen_id))

            child.genotype.chromosome_v.gens[mut_gen_id].point['z'] = z_drawn
        else:
            """ Children is not about to be mutated in vertical, just pass """
            pass

        """ Tangent angle h mutation """
        mutate_mask = np.random.choice(10, 1)[0]
        """ Drawning 0 is 10% chance - mutations probability """
        if mutate_mask == 0:

            """ Draw new angle value and replace """
            child.genotype.chromosome_h.init_tangent_rand()

            print('\t\t\tMutation will be applied for ANGLE H: {}'
                  .format(child.genotype.chromosome_h.gen_init_tangent))

        else:
            """ Children is not about to be mutated in angle h, just pass """
            pass

        """ Tangent angle v mutation """
        mutate_mask = np.random.choice(10, 1)[0]
        """ Drawning 0 is 10% chance - mutations probability """
        if mutate_mask == 0:

            """ Draw new angle value and replace """
            child.genotype.chromosome_v.init_tangent_rand()

            print('\t\t\tMutation will be applied for ANGLE V: {}'
                  .format(child.genotype.chromosome_v.gen_init_tangent))

        else:
            """ Children is not about to be mutated in angle v, just pass """
            pass

        return child


def is_point_valid(x: int, y: int) -> bool:
    """
    Check if drawn coordinates are in map range and if corrensponding z value in not np.inf
    Params:                                                                     type:
    :param x: Drawn x coordinate                                                int
    :param y: Drawn y coordinate                                                int
    :return: True if drawn point is in map range, else False                    bool
    """
    if x in X_DATA_INV and y in Y_DATA_INV:
        return False
    else:
        return True


def sort_nearest_order(xvals: List[Union[int, float]], yvals: List[Union[int, float]], end_point: Tuple[int, int]) ->\
        (List[Tuple[Union[int, float]]]):
    distances = []
    for point_id in range(len(xvals)):
        test_point = sp.Point2D(xvals[point_id], yvals[point_id])
        x_end, y_end = end_point
        end_point = sp.Point2D(x_end, y_end)
        distances.append(float(test_point.distance(end_point)))

    return [(xvals[point_id], yvals[point_id]) for distance, point_id in
            sorted(zip(distances, range(len(xvals))), reverse=True)]


def print_population_info(title: str, pop: np.array) -> None:
    pop = deepcopy(pop)
    print('DB:', title)
    gens_descs_h = [[(math.floor(gen.point['x']), math.floor(gen.point['y'])) for
                     gen in ind.genotype.chromosome_h.gens] for ind in pop]
    gens_descs_v = [[(math.floor(gen.point['d']), math.floor(gen.point['z'])) for
                     gen in ind.genotype.chromosome_v.gens] for ind in pop]

    crcs_h = [ind.genotype.checksum_h for ind in pop]
    crcs_v = [ind.genotype.checksum_v for ind in pop]

    tangs_h = [round(float(ind.genotype.chromosome_h.gen_init_tangent), 2) for ind in pop]
    tangs_v = [round(float(ind.genotype.chromosome_v.gen_init_tangent), 2) for ind in pop]

    i = 0
    for index in range(len(pop)):
        print('\t{}\tH:'.format(i), gens_descs_h[index], '[{:4.2f}]'.format(tangs_h[index]), 'CRC:', crcs_h[index])
        print('\t\tV:', gens_descs_v[index], '[{:4.2f}]'.format(tangs_v[index]), 'CRC:', crcs_v[index])
        print('\t\tF:', format_fitness(fitness=pop[index].fenotype.fitness_val))
        i += 1


def format_fitness(fitness: float) -> str:
    f = str(math.floor(fitness))
    l = len(f)
    n = math.floor(l / 3)

    if n > 0 and l > 3:
        r = l % 3
        if r == 0:
            n_dots = n - 1
        else:
            n_dots = n

        for i in range(n_dots):
            f = f[:-3 * (i + 1) - i] + "'" + f[-3 * (i + 1) - i:]

    return f




class GAModel:
    def __init__(self):
        global X_DATA_FIL, Y_DATA_FIL, Z_DATA_FIL
        X_DATA_FIL, Y_DATA_FIL, Z_DATA_FIL = read_filtered_data()
        global X_DATA_INV, Y_DATA_INV
        X_DATA_INV, Y_DATA_INV = read_invalid_data()

        self.population = Population(pop_size=POPULATION_SIZE)
        self.population.initialize_random()
        best_ind = self.population.get_best_individual()
        print('\t^ Best fitness:', format_fitness(fitness=best_ind.fenotype.fitness_val), '\n')

    def evaluate(self):
        best_ind = None

        fitness_hist = []

        for i in range(GENERATIONS_NUM):
            print('Create new generation:', i)
            self.population.create_new_generation()

            best_ind = deepcopy(self.population.get_best_individual())
            fitness = best_ind.fenotype.fitness_val
            print('\t^ Best fitness:', format_fitness(fitness=fitness), '\n')
            fitness_hist.append(fitness)

            plt.plot(range(len(fitness_hist)), fitness_hist)
            plt.title('fitness: i')
            plt.grid()
            plt.xlabel('Numer iteracji')
            plt.ylabel('Koszt [euro]')
            # plt.show()
            plt.savefig('output\\' + '_fitness_' + str(i) + '.png')
            plt.close()

        """ Plot best route """
        best_route_desc_h = best_ind.fenotype.get_route_desc(plane=Plane.HORIZONTAL)
        best_route_len_h = best_ind.fenotype.route_len_h
        best_p_descs_h = best_ind.genotype.get_points_descs(plane=Plane.HORIZONTAL)
        plot_route_2d(plane=Plane.HORIZONTAL, route_desc=best_route_desc_h, route_len=best_route_len_h,
                      p_dicts=best_p_descs_h, title='h: end best')

        dicretized_points = best_ind.dicretize_route()
        landform_hights = best_ind.get_landform_route_heights(route_points=dicretized_points)
        d_points_disc = [elem['d'] for elem in dicretized_points]
        landform = (d_points_disc, landform_hights)

        best_route_desc_v = best_ind.fenotype.get_route_desc(plane=Plane.VERTICAL)
        best_route_len_v = best_ind.fenotype.route_len_total
        plot_route_2d(plane=Plane.VERTICAL, route_desc=best_route_desc_v, route_len=best_route_len_v,
                      p_dicts=best_ind.map_genotype_v_to_p_dicts(), landform=landform, title='v: end best',
                      z_min=min(landform_hights), z_max=max(landform_hights))

        plt.plot(range(len(fitness_hist)), fitness_hist)
        plt.title('fitness: end')
        plt.grid()
        plt.xlabel('Numer iteracji')
        plt.ylabel('Koszt [euro]')
        # plt.show()
        plt.savefig('output\\' + '_fitness_change' + '.png')
        plt.close()
