from HypTypes import *
from copy import copy, deepcopy
from typing import Dict, Tuple, List, Union
from HypGeo import get_route_description, plot_route_2d
import collections
from HypData import read_invalid_data
import sympy as sp
import math


class Gen:
    """ Definition of Gen - coordinates of single point """
    def __init__(self, plane: Plane) -> None:
        self.point = deepcopy(POINT_DEF[plane])


class Chromosome:
    """ Definition of Chromosome - collection of Gens """
    def __init__(self, plane: Plane) -> None:
        self.gens = np.array([Gen(plane=plane) for i in range(CHROMOSOME_SIZE)])
        self.gen_init_tangent = np.inf

    def rand_init_tangent(self, plane: Plane):
        if plane==Plane.HORIZONTAL:
            """ Initialize tangent in first point (for horizontal drawn from range (0, pi/2))"""
            self.gen_init_tangent = sp.pi / 2 * np.random.sample()
        else:
            raise NotImplementedError


class Genotype:
    """ Definition of Genotype - collection of Chromosomes """
    def __init__(self) -> None:
        self.chromosome_h = Chromosome(plane=Plane.HORIZONTAL)
        # self.chromosome_v = Chromosome(plane=Plane.VERTICAL)
        self.checksum_h = np.inf

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

        """ Sort points to be from most far away to the closest to end point """
        values_x = deepcopy([gen.point['x'] for gen in gens[1:-1]])
        values_y = deepcopy([gen.point['y'] for gen in gens[1:-1]])
        values_x.sort()
        values_y.sort()

        # sorted_nearest = sort_nearest_order(xvals=values_x, yvals=values_y,
        #                                     end_point=(END_POINT['x'], END_POINT['y']))

        # for point in sorted_nearest:
        #     (gens[sorted_nearest.index(point) + 1]).point['x'], (gens[sorted_nearest.index(point) + 1]).point['y'] =\
        #         point

        for gen_id in range(1, len(gens) - 1):
            (gens[gen_id]).point['x'] = values_x[gen_id - 1]
            (gens[gen_id]).point['y'] = values_y[gen_id - 1]

        # for gen_id in range(1, len(gens) - 1):
        #     x_sorted, y_sorted = sorted_nearest[gen_id - 1]
        #     (gens[gen_id]).point['x'] = x_sorted
        #     (gens[gen_id]).point['y'] = y_sorted

        """ Initialize tangent in first point (for horizontal drawn from range (0, pi/2))"""
        self.chromosome_h.gen_init_tangent = sp.pi/2 * np.random.sample()

        """ Calculate checksum - only intermediate points plus tangent value times 100 """
        self.checksum_h = math.floor(sum(values_x + values_y)) + float(self.chromosome_h.gen_init_tangent) * 100

        # print('Genotype: randomly initialized horizontal: {}, {}'.
        #       format([gen.point for gen in gens], self.chromosome_h.gen_init_tangent))

    # def init_vertical(self, route_len_h: float) -> None:
    #     """
    #     Initialize Gens (Points) for vertical movement
    #     Params:                                                                     type:
    #     :param route_len_h: route length in horizontal plane                        float
    #     :return: None
    #     """
    #
    #     """ Initialize start point and end point coordinates """
    #     gens = self.chromosome_v.gens
    #     gens[0].point['z'] = START_POINT['z']
    #     gens[0].point['d'] = 0
    #     gens[-1].point['z'] = END_POINT['z']
    #     gens[-1].point['d'] = route_len_h
    #
    #     """ Draw intermediate points coordinates """
    #     z_drawns = np.random.choice([i for i in range(ROUTE_MIN_HIGHT, ROUTE_MAX_HIGHT)], CHROMOSOME_SIZE - 2)
    #     d_drawns = np.random.choice([i for i in range(1, math.floor(route_len_h))], CHROMOSOME_SIZE - 2)
    #
    #     """ Check if d values are unique """
    #     seen = set()
    #     if any(i in seen or seen.add(i) for i in d_drawns):
    #         raise ValueError('Values in d_drawns not unique: {}!'.format(d_drawns))
    #
    #     """ Sort values by growing distance d """
    #     df_vert = pd.DataFrame({'z': z_drawns, 'd': d_drawns})
    #     df_vert = df_vert.sort_values(by=['d'], ignore_index=True)
    #
    #     """ Set intermediate points """
    #     for gen_id in range(1, CHROMOSOME_SIZE - 1):
    #         gens[gen_id].point['z'] = df_vert['z'][gen_id - 1]
    #         gens[gen_id].point['d'] = df_vert['d'][gen_id - 1]
    #
    #     """ Initialize tangent in first point (for vertical drawn from range (-pi/8, pi/8))"""
    #     self.chromosome_v.gen_init_tangent = sp.pi/4 * np.random.sample() - sp.pi/8
    #
    #     print('Genotype: randomly initialized vertical: {}'.format([gen.point for gen in gens]))

    def get_points_descs(self, plane: Plane) -> List[Dict[str, int]]:
        if plane == Plane.HORIZONTAL:
            return [{'x': gen.point['x'], 'y': gen.point['y']} for gen in self.chromosome_h.gens]
        # elif plane == Plane.VERTICAL:
        #     return [{'z': gen.point['z'], 'd': gen.point['d']} for gen in self.chromosome_v.gens]

    def init_checksum_h(self):
        gens = self.chromosome_h.gens
        values_x = [gen.point['x'] for gen in gens[1:-1]]
        values_y = [gen.point['y'] for gen in gens[1:-1]]
        self.checksum_h = math.floor(sum(values_x + values_y)) + float(self.chromosome_h.gen_init_tangent) * 100
        return copy(self.checksum_h)


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
        # self.route_vertical = RouteDesc()
        self.route_horizontal_length = np.inf
        # self.route_total_length = np.inf
        self.route_horizontal_desc = np.inf

    def init_horizontal(self, p_dicts: List[Dict[str, int]], init_tangent: float) -> None:
        route_desc, route_len_h = get_route_description(plane=Plane.HORIZONTAL, p_dicts=p_dicts,
                                                        init_tangent=init_tangent)
        arcs = self.route_horizontal.arcs
        for (arc, arc_desc) in zip(arcs, route_desc):
            arc.point_cR = arc_desc['point_cR']
            arc.circle_radius_len = arc_desc['circle_radius_len']
            arc.ray_RA = arc_desc['ray_RA']
            arc.arc_rad_len = arc_desc['arc_rad_len']
            arc.ray_tangent_B = arc_desc['ray_tangent_B']
            arc.segment = arc_desc['segment']

        self.route_horizontal_length = route_len_h
        self.route_horizontal_desc = route_desc
        # print('Fenotype: initialized horizontal: {}'.format(self.route_horizontal_desc))

        if PLOT_INIT:
            plot_route_2d(plane=Plane.HORIZONTAL, route_desc=route_desc, route_len=route_len_h, p_dicts=p_dicts)

    # def init_vertical(self, p_dicts: List[Dict[str, int]], init_tangent: float) -> None:
    #     route_desc, route_len = get_route_description(plane=Plane.VERTICAL, p_dicts=p_dicts, init_tangent=init_tangent)
    #     arcs = self.route_vertical.arcs
    #     for (arc, arc_desc) in zip(arcs, route_desc):
    #         arc.point_cR = arc_desc['point_cR']
    #         arc.ray_RA = arc_desc['ray_RA']
    #         arc.arc_rad_len = arc_desc['arc_rad_len']
    #         arc.ray_tangent_B = arc_desc['ray_tangent_B']
    #         arc.segment = arc_desc['segment']
    #
    #     self.route_total_length = route_len
    #     print('Fenotype: initialized vertical: {}'.format(self.route_desc))
    #
    #     if PLOT_INIT:
    #         plot_route_2d(plane=Plane.VERTICAL, route_desc=route_desc, route_len=self.route_horizontal_length,
    #                       p_dicts=p_dicts)

    def get_route_desc(self):
        return self.route_horizontal_desc

    def fitness(self) -> float:
        return self.route_horizontal_length


class Individual:
    def __init__(self):
        self.genotype = Genotype()
        self.fenotype = Fenotype()

    def initialize_random(self, invalid_coordinates: Tuple[List[int], List[int]]):
        # print('Individual: started random initialization')
        self.genotype.init_horizontal(invalid_coordinates=invalid_coordinates)
        self.fenotype.init_horizontal(p_dicts=self.genotype.get_points_descs(plane=Plane.HORIZONTAL),
                                      init_tangent=self.genotype.chromosome_h.gen_init_tangent)

        # self.genotype.init_vertical(route_len_h=self.fenotype.route_horizontal_length)
        # self.fenotype.init_vertical(p_dicts=self.genotype.get_points_descs(plane=Plane.VERTICAL),
        #                             init_tangent=self.genotype.chromosome_v.gen_init_tangent)
        # print('Individual: ended random initialization')

    def calculate_cost(self) -> float:
        return self.fenotype.route_horizontal_length


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

    def create_new_generation(self) -> None:
        #
        print('\n', '********************')
        print('DB:', 'start gen new pop')
        old_pop = [[(math.floor(gen.point['x']), math.floor(gen.point['y'])) for gen in ind.genotype.chromosome_h.gens] for ind in self.individuals]
        print('DB:', 'old pop:')
        i = 0
        for ind in old_pop:
            print('\t', i, ind, 'CRC:', math.floor(sum([x + y for (x, y) in ind]) / 1000))
            i += 1
        #

        """ Create new offspring """
        offspring = np.array([Individual() for i in range(len(self.individuals))])

        #
        i = 0
        print('DB:', 'ofspring part 1')
        offsp1 = [[(math.floor(gen.point['x']) if gen.point['x'] != np.inf else 999999, math.floor(gen.point['y']) if gen.point['x'] != np.inf else 999999) for gen in ind.genotype.chromosome_h.gens] for ind in offspring]
        for ind in offsp1:
            print('\t', i, ind, 'CRC:', math.floor(sum([x + y for (x, y) in ind]) / 1000))
            i += 1
        #

        """ select 30% best parents """
        best_to_worst_parents = self.__get_best_individuals_sorted()
        # print('***', [p.calculate_cost() for p in best_to_worst_parents])

        #
        i = 0
        print('DB:', 'best parents')
        best_par = [[(math.floor(gen.point['x']) if gen.point['x'] != np.inf else 999999, math.floor(gen.point['y']) if gen.point['x'] != np.inf else 999999) for gen in ind.genotype.chromosome_h.gens] for ind in best_to_worst_parents]
        for ind in best_par:
            print('\t', i, ind, 'CRC:', math.floor(sum([x + y for (x, y) in ind]) / 1000))
            i += 1
        #

        offspring[NEWPOP_BEST_PARENTS_START:NEWPOP_BEST_PARENTS_NUM] =\
            deepcopy(best_to_worst_parents[0:NEWPOP_BEST_PARENTS_NUM])
        # print('Population: selected best parents')

        #
        i = 0
        print('DB:', 'ofspring part 2')
        offsp2 = [[(math.floor(gen.point['x']) if gen.point['x'] != np.inf else 999999, math.floor(gen.point['y']) if gen.point['x'] != np.inf else 999999) for gen in ind.genotype.chromosome_h.gens] for ind in offspring]
        for ind in offsp2:
            print('\t', i, ind, 'CRC:', math.floor(sum([x + y for (x, y) in ind]) / 1000))
            i += 1
        #

        """ use best parents to cross over offspring and generate 50% of them """
        children = self.__cross_over2(num_of_ind=NEWPOP_CHILDREN_NUM, parents_sorted=best_to_worst_parents)
        for individual in children:
            individual.fenotype.init_horizontal(p_dicts=individual.genotype.get_points_descs(plane=Plane.HORIZONTAL),
                                                init_tangent=individual.genotype.chromosome_h.gen_init_tangent)

        offspring[NEWPOP_CHILDREN_START:NEWPOP_CHILDREN_START + NEWPOP_CHILDREN_NUM] = deepcopy(children)
        # print('Population: produced new children')

        #
        i = 0
        print('DB:', 'ofspring part 3')
        offsp3 = [[(math.floor(gen.point['x']) if gen.point['x'] != np.inf else 999999,
                    math.floor(gen.point['y']) if gen.point['x'] != np.inf else 999999) for gen in
                   ind.genotype.chromosome_h.gens] for ind in offspring]
        for ind in offsp3:
            print('\t', i, ind, 'CRC:', math.floor(sum([x + y for (x, y) in ind]) / 1000))
            i += 1
        #

        """ create random 20% individuals """
        random_pop = np.array([Individual() for i in range(NEWPOP_RANDOM_NUM)])
        x_inv, y_inv = read_invalid_data()
        for individual in random_pop:
            individual.initialize_random(invalid_coordinates=(x_inv, y_inv))

        #
        i = 0

        print('DB:', 'randoms')
        randoms = [[(math.floor(gen.point['x']) if gen.point['x'] != np.inf else 999999, math.floor(gen.point['y']) if gen.point['x'] != np.inf else 999999) for gen in ind.genotype.chromosome_h.gens] for ind in random_pop]
        for ind in randoms:
            crc = math.floor(sum([x + y for (x, y) in ind]) / 1000)
            print('\t', i, ind, 'CRC:', crc)
            i += 1
        #
        offspring[NEWPOP_RANDOM_START:NEWPOP_RANDOM_START + NEWPOP_RANDOM_NUM] = deepcopy(random_pop)
        # print('Population: produced randoms, len: {}'.format(len(random_pop)))

        #
        i = 0
        crcs = []
        print('DB:', 'ofspring part 4')
        offsp4 = [[(math.floor(gen.point['x']) if gen.point['x'] != np.inf else 999999,
                    math.floor(gen.point['y']) if gen.point['x'] != np.inf else 999999) for gen in
                   ind.genotype.chromosome_h.gens] for ind in offspring]
        for ind in offsp4:
            crc = math.floor(sum([x + y for (x, y) in ind]) / 1000)
            print('\t', i, ind, 'CRC:', crc)
            crcs.append(crc)
            i += 1
        print('Unique len:', len([item for item, count in collections.Counter(crcs).items() if count == 1]))
        if len([item for item, count in collections.Counter(crcs).items() if count == 1]) != 10:
            print('DUPLICATES !!! !!! !!!')
        #

        """ Set new population """
        for (individual, index) in zip(self.individuals, range(len(self.individuals))):
            individual.genotype = deepcopy(offspring[index].genotype)
            individual.fenotype = deepcopy(offspring[index].fenotype)
        # print('Population: made offspring new population')
        print('####################', '\n')

    def __get_best_individuals_sorted(self) -> np.array:
        individuals = deepcopy(self.individuals)
        costs_desc = deepcopy([{'index': index, 'cost': individual.calculate_cost()} for
                              (index, individual) in zip(range(len(individuals)), individuals)])
        # print([item['cost'] for item in costs_desc])

        sorted_indexes = [index for (cost, index) in sorted(zip([item['cost'] for item in costs_desc],
                                                                [item['index'] for item in costs_desc]))]

        # self.individuals = np.array([self.individuals[index] for index in sorted_indexes])
        # print(sorted([item['cost'] for item in costs_desc]))
        return np.array([individuals[index] for index in sorted_indexes])

    def get_best_individual(self) -> Individual:
        return (self.__get_best_individuals_sorted())[0]

    def __cross_over2(self, num_of_ind: int, parents_sorted: np.array) -> np.array:
        """ """

        """ Create new empty population - children """
        children = np.array([Individual() for i in range(num_of_ind)])
        parents = deepcopy(parents_sorted)
        child_id = 0

        """ Create list of mating ranges for every parent - every parent index from list parents indicates
            top value of selection range, bottom selection range is previous value in list """
        parents_mating_points = []
        for parent_id in range(len(parents)):
            parents_mating_points.append(math.ceil(parents[0].fenotype.fitness() / parents[parent_id].fenotype.fitness() * 100)
                                         + parents_mating_points[parent_id -1] if parent_id > 0 else 0)
        mating_points_sum = parents_mating_points[-1]

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

            """ Check if chosen parents are not the same Individual """
            if parent1.genotype.checksum_h != parent2.genotype.checksum_h:
                """ Generate crossing mask and make sure it is not uniform """
                while True:
                    cross_mask = np.random.choice([0, 1], CHROMOSOME_SIZE - 2)
                    if not sum(cross_mask) in [0, len(cross_mask)]:
                        break

                """ Generate mask for choosing init tangent in starting point"""
                tang_mask = np.random.choice([0, 1], 1)[0]

                """ Start and end point is always the same """
                child.genotype.chromosome_h.gens[0] = deepcopy(parent1.genotype.chromosome_h.gens[0])
                child.genotype.chromosome_h.gens[-1] = deepcopy(parent1.genotype.chromosome_h.gens[-1])

                """ Copy intermediate gens according to crossing mask """
                for i in range(1, CHROMOSOME_SIZE - 1):
                    child.genotype.chromosome_h.gens[i] = \
                        deepcopy((parent1.genotype.chromosome_h.gens[i]) if (cross_mask[i - 1] == 0)
                                 else (parent2.genotype.chromosome_h.gens[i]))

                child.genotype.chromosome_h.gen_init_tangent = deepcopy(
                    parent1.genotype.chromosome_h.gen_init_tangent if (tang_mask == 0)
                    else parent2.genotype.chromosome_h.gen_init_tangent)

                """ Check if created individual does not already exist """
                crc_child = child.genotype.init_checksum_h()
                if crc_child not in [ind.genotype.checksum_h for ind in parents]:
                    """ Add individual to new population """
                    children[child_id] = deepcopy(child)
                    child_id += 1

                    if child_id == num_of_ind:
                        """ Produced all children, return """
                        return children

                    else:
                        """ Not all children produced yet - go to beginning of the loop """
                        pass
                else:
                    """ Individual already exists in population - go to beginning of the loop """
                    pass
            else:
                """ Invalid parents drawing - go to beginning of the loop """
                pass









    def __cross_over(self, num_of_ind: int, parents_sorted: np.array) -> np.array:
        children = np.array([Individual() for i in range(num_of_ind)])
        parents = deepcopy(parents_sorted)

        #
        print('DB:', 'crossed children')
        db_index = 0
        #

        """ cross horizontal """
        selector_id = 0
        for individual in children:
            """ make sure mask elements are different so it will not just copy parent """
            while True:
                cross_mask = np.random.choice([0, 1], CHROMOSOME_SIZE - 2)
                mask_begin = cross_mask[0]
                mask_uniform = True
                for mask_elem in cross_mask[1:]:
                    if mask_elem != mask_begin:
                        mask_uniform = False
                        break
                if mask_uniform is False:
                    break
            # print(cross_mask)

            parent1_id, parent2_id = CROSS_SELECT5[selector_id]
            parent1, parent2 = parents[parent1_id], parents[parent2_id]
            individual.genotype.chromosome_h.gens[0] = parent1.genotype.chromosome_h.gens[0]
            individual.genotype.chromosome_h.gens[-1] = parent1.genotype.chromosome_h.gens[-1]
            for i in range(1, CHROMOSOME_SIZE - 1):
                individual.genotype.chromosome_h.gens[i] = deepcopy((parent1.genotype.chromosome_h.gens[i])
                                                                    if (cross_mask[i - 1] == 0)
                                                                    else (parent2.genotype.chromosome_h.gens[i]))
            tangent_mask = np.random.choice([0, 1], 1)[0]

            #
            single_ch = [(math.floor(gen.point['x']) if gen.point['x'] != np.inf else 999999, math.floor(gen.point['y']) if gen.point['x'] != np.inf else 999999) for gen in individual.genotype.chromosome_h.gens]
            par1 = [(math.floor(gen.point['x']) if gen.point['x'] != np.inf else 999999, math.floor(gen.point['y']) if gen.point['x'] != np.inf else 999999) for gen in parent1.genotype.chromosome_h.gens]
            par2 = [(math.floor(gen.point['x']) if gen.point['x'] != np.inf else 999999, math.floor(gen.point['y']) if gen.point['x'] != np.inf else 999999) for gen in parent2.genotype.chromosome_h.gens]
            print('\t', db_index, single_ch, 'CRC:', math.floor(sum([x + y for (x, y) in single_ch]) / 1000))
            print('\t', 'mask:', cross_mask, [tangent_mask], 'parents:', parent1_id, parent2_id)
            print('\t\t', '0:', par1, 'CRC:', math.floor(sum([x + y for (x, y) in par1]) / 1000))
            print('\t\t', '1:', par2, 'CRC:', math.floor(sum([x + y for (x, y) in par2]) / 1000))
            db_index += 1
            #


            # print(tangent_mask)
            # print(parent1_id, parent2_id)
            individual.genotype.chromosome_h.gen_init_tangent =\
                deepcopy((parent1.genotype.chromosome_h.gen_init_tangent)
                         if (tangent_mask == 0)
                         else (parent2.genotype.chromosome_h.gen_init_tangent))

            """ cross vertical """
            # cross_mask = np.random.choice([0, 1], CHROMOSOME_SIZE)
            # for i in range(CHROMOSOME_SIZE):
            #     individual.genotype.chromosome_v.gens[i] = parent1.genotype.chromosome_v.gens[i] if cross_mask[i] == 0 \
            #         else parent2.genotype.chromosome_v.gens[i]
            # TODO get tangent from one parent

            selector_id += 1

        return children


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
