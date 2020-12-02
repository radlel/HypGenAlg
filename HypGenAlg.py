from HypTypes import *
from copy import copy, deepcopy
from typing import Dict, Tuple, List, Union
from HypGeo import get_route_description, plot_route_2d
from HypData import read_invalid_data
import sympy as sp
import pandas as pd
import math


""" Genetic Algorithm parameters """
CHROMOSOME_SIZE = 5
GENERATIONS_NUM = 40

NEWPOP_BEST_PARENTS_NUM = 15
NEWPOP_CHILDREN_NUM = 25
NEWPOP_RANDOM_NUM = 10

POPULATION_SIZE = NEWPOP_BEST_PARENTS_NUM + NEWPOP_CHILDREN_NUM + NEWPOP_RANDOM_NUM

NEWPOP_BEST_PARENTS_START = 0
NEWPOP_CHILDREN_START = NEWPOP_BEST_PARENTS_NUM
NEWPOP_RANDOM_START = NEWPOP_BEST_PARENTS_NUM + NEWPOP_CHILDREN_NUM

ROUTE_RESOLUTION = 1000
MATING_POINTS_MAX = 100


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
            """ Initialize tangent in first point (for horizontal drawn from range (0, pi/2)) """
            self.gen_init_tangent = sp.pi/2 * np.random.sample()
        elif self.plane == Plane.VERTICAL:
            """ Initialize tangent in first point (for horizontal drawn from range (-pi/8, pi/8)) """
            self.gen_init_tangent = sp.pi/4 * np.random.sample() - sp.pi/8
        else:
            raise ValueError('Invalid plane parameter value! {}'.format(self.plane))


class Genotype:
    """ Definition of Genotype - collection of Chromosomes """
    def __init__(self) -> None:
        self.chromosome_h = Chromosome(plane=Plane.HORIZONTAL)
        self.checksum_h = np.inf
        self.chromosome_v = Chromosome(plane=Plane.VERTICAL)
        self.checksum_v = np.inf

    def init_horizontal(self, invalid_coordinates: Tuple[List[int], List[int]]) -> None:
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
                x_drawn = np.random.choice([i for i in range(MAP_LIMIT['xmin'] + DIST10KM,
                                                             MAP_LIMIT['xmax'] - DIST10KM)], 1)[0]
                y_drawn = np.random.choice([i for i in range(MAP_LIMIT['ymin'] + DIST10KM,
                                                             MAP_LIMIT['ymax'] - DIST10KM)], 1)[0]

                """ Check if selected point has proper z value (not np.inf) """
                if is_point_valid(x_drawn=x_drawn, y_drawn=y_drawn, invalid_coordinates=invalid_coordinates):
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

    def init_vertical(self) -> None:
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
        d_drawns = np.random.choice(ROUTE_RESOLUTION, CHROMOSOME_SIZE - 2)

        """ Check if d values are unique """
        # TODO - check if d values are unique

        """ Sort values by growing distance d """
        df_vert = pd.DataFrame({'z': z_drawns, 'd': d_drawns})
        df_vert = deepcopy(df_vert.sort_values(by=['d'], ignore_index=True))

        """ Set intermediate points """
        for gen_id in range(1, CHROMOSOME_SIZE - 1):
            gens[gen_id].point['z'] = df_vert['z'][gen_id - 1]
            gens[gen_id].point['d'] = df_vert['d'][gen_id - 1]

        """ Initialize tangent in first point (for vertical drawn from range (-pi/8, pi/8))"""
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
            plot_route_2d(plane=Plane.HORIZONTAL, route_desc=route_desc, route_len=route_len, p_dicts=p_dicts)

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
            plot_route_2d(plane=Plane.VERTICAL, route_desc=route_desc, route_len=route_len, p_dicts=p_dicts)

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

    def initialize_random(self, invalid_coordinates: Tuple[List[int], List[int]]) -> None:
        self.genotype.init_horizontal(invalid_coordinates=invalid_coordinates)
        self.fenotype.init_horizontal(p_dicts=self.genotype.get_points_descs(plane=Plane.HORIZONTAL),
                                      init_tangent=self.genotype.chromosome_h.gen_init_tangent)

        self.genotype.init_vertical()
        self.fenotype.init_vertical(p_dicts=self.map_genotype_v_to_p_dicts(),
                                    init_tangent=self.genotype.chromosome_v.gen_init_tangent)

    def map_genotype_v_to_p_dicts(self) -> List[Dict[str, int]]:
        points_v = deepcopy(self.genotype.get_points_descs(plane=Plane.VERTICAL))
        route_len_h = copy(self.fenotype.route_len_h)
        return [{'z': point['z'], 'd': math.floor(point['d'] / ROUTE_RESOLUTION * route_len_h)} for point in points_v]

    def fitness(self) -> float:
        # TODO - route len as a cost only temporarly !!!
        return self.fenotype.route_len_total


class Population:
    """ Definition of Population - collection of Genotypes """
    def __init__(self, pop_size=POPULATION_SIZE) -> None:
        self.individuals = np.array([Individual() for i in range(pop_size)])

    def initialize_random(self) -> None:
        """
        Initialize Individuals with random values
        Params:                                                                     type:
        :return: None                                                               None
        """
        print('Population: started random initialization')
        """ Get invalid coordinates so Genotype can check if drawn point is correct """
        x_inv, y_inv = read_invalid_data()

        for individual in self.individuals:
            print('Population: initialize Indywidual {}'.format((list(self.individuals)).index(individual)))
            individual.initialize_random(invalid_coordinates=(x_inv, y_inv))
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
        children = self.__cross_over(parents_sorted=best_to_worst_parents)
        for individual in children:
            individual.fenotype.init_horizontal(p_dicts=individual.genotype.get_points_descs(plane=Plane.HORIZONTAL),
                                                init_tangent=individual.genotype.chromosome_h.gen_init_tangent)

            individual.fenotype.init_vertical(p_dicts=individual.map_genotype_v_to_p_dicts(),
                                              init_tangent=individual.genotype.chromosome_v.gen_init_tangent)

        offspring[NEWPOP_CHILDREN_START:NEWPOP_CHILDREN_START + NEWPOP_CHILDREN_NUM] = deepcopy(children)

        print_population_info(title='Children after crossing over', pop=children)

        """ create random 20% individuals """
        random_pop = Population(pop_size=NEWPOP_RANDOM_NUM)
        random_pop.initialize_random()
        offspring[NEWPOP_RANDOM_START:NEWPOP_RANDOM_START + NEWPOP_RANDOM_NUM] = deepcopy(random_pop.individuals)

        """ Set new population """
        for (individual, index) in zip(self.individuals, range(len(self.individuals))):
            individual.genotype = deepcopy(offspring[index].genotype)
            individual.fenotype = deepcopy(offspring[index].fenotype)

        print_population_info(title='New population', pop=self.individuals)

        print('...FINISHED CREATION OF NEW GENERATION')

    def __get_best_individuals_sorted(self) -> np.array:
        individuals = deepcopy(self.individuals)
        costs_desc = deepcopy([{'index': index, 'cost': individual.fitness()} for
                              (index, individual) in zip(range(len(individuals)), individuals)])
        # print([item['cost'] for item in costs_desc])

        sorted_indexes = [index for (cost, index) in sorted(zip([item['cost'] for item in costs_desc],
                                                                [item['index'] for item in costs_desc]))]

        # self.individuals = np.array([self.individuals[index] for index in sorted_indexes])
        # print(sorted([item['cost'] for item in costs_desc]))
        return np.array([individuals[index] for index in sorted_indexes])

    def get_best_individual(self) -> Individual:
        return (self.__get_best_individuals_sorted())[0]

    def __cross_over(self, parents_sorted: np.array) -> np.array:
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
            parents_mating_points.append(math.ceil(parents[0].fitness() /
                                                   parents[parent_id].fitness() *
                                                   MATING_POINTS_MAX)
                                         + (parents_mating_points[parent_id -1] if parent_id > 0 else 0))
        mating_points_sum = parents_mating_points[-1]
        print('\t\t\t\tParents mating points:', parents_mating_points)

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

                child.genotype.chromosome_h.gen_init_tangent = deepcopy(
                    parent1.genotype.chromosome_h.gen_init_tangent if (tang_mask_h == 0)
                    else parent2.genotype.chromosome_h.gen_init_tangent)

                child.genotype.chromosome_v.gen_init_tangent = deepcopy(
                    parent1.genotype.chromosome_v.gen_init_tangent if (tang_mask_v == 0)
                    else parent2.genotype.chromosome_v.gen_init_tangent)

                """ Check if created individual does not already exist """
                crc_h_child = child.genotype.init_checksum(plane=Plane.HORIZONTAL)
                crc_v_child = child.genotype.init_checksum(plane=Plane.VERTICAL)
                if ((crc_h_child not in [ind.genotype.checksum_h for ind in parents] and
                     (crc_v_child not in [ind.genotype.checksum_v for ind in parents]))):
                    """ Add individual to new population """
                    children[child_id] = deepcopy(child)
                    print('\t\t\t\tCrossed over child:', child_id)
                    child_id += 1

                    if child_id == NEWPOP_CHILDREN_NUM:
                        """ Produced all children, return """
                        print('\t\t\t...Ended cross over')
                        return children

                    else:
                        """ Not all children produced yet - go to beginning of the loop """
                        pass
                else:
                    """ Individual already exists in population - go to beginning of the loop """
                    print('\t\t\t*** Individual with at least one the same CRC already exists in population, draw '
                          'another parents pair')
            else:
                """ Invalid parents drawing - go to beginning of the loop """
                print('\t\t\t*** Parents has at least one same CRC, draw another parents pair')


def is_point_valid(x_drawn: int, y_drawn: int, invalid_coordinates: Tuple[List[int], List[int]]) -> bool:
    """
    Check if drawn coordinates are in map range and if corrensponding z value in not np.inf
    Params:                                                                     type:
    :param x_drawn: Drawn x coordinate                                          int
    :param y_drawn: Drawn y coordinate                                          int
    :param invalid_coordinates:
    :return: True if drawn point is in map range, else False                    bool
    """
    x_inv, y_inv = invalid_coordinates
    if x_drawn in x_inv and y_drawn in y_inv:
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
        print('\t\tF:', math.floor(pop[index].fitness()) / 1000)
        i += 1


class GAModel:
    def __init__(self):
        self.population = Population(pop_size=POPULATION_SIZE)
        self.population.initialize_random()
        best_ind = self.population.get_best_individual()
        print('\t^ Best fitness:', math.floor(best_ind.fitness()) / 1000, '\n')

    def evaluate(self):
        best_ind = None
        for i in range(GENERATIONS_NUM):
            print('Create new generation:', i)
            self.population.create_new_generation()

            best_ind = deepcopy(self.population.get_best_individual())
            print('\t^ Best fitness:', math.floor(best_ind.fitness()) / 1000, '\n')

        """ Plot best route """
        best_route_desc_h = best_ind.fenotype.get_route_desc(plane=Plane.HORIZONTAL)
        best_route_len_h = best_ind.fenotype.route_len_h
        best_p_descs_h = best_ind.genotype.get_points_descs(plane=Plane.HORIZONTAL)
        plot_route_2d(plane=Plane.HORIZONTAL, route_desc=best_route_desc_h, route_len=best_route_len_h,
                      p_dicts=best_p_descs_h)

        best_route_desc_v = best_ind.fenotype.get_route_desc(plane=Plane.VERTICAL)
        best_route_len_v = best_ind.fenotype.route_len_total
        plot_route_2d(plane=Plane.VERTICAL, route_desc=best_route_desc_v, route_len=best_route_len_v,
                      p_dicts=best_ind.map_genotype_v_to_p_dicts())
