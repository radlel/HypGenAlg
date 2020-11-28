from HypTypes import *
from copy import deepcopy
from typing import Dict
import pandas as pd
from HypGeo import evaluate_route_len, plot_route_2d
import math
from HypData import read_filtered_data


class Gen:
    """ Definition of Gen - coordinates of single point """
    def __init__(self, plane: Plane) -> None:
        self.point = deepcopy(POINT_DEF[plane])


class Chromosome:
    """ Definition of Chromosome - collection of Gens """
    def __init__(self, plane: Plane) -> None:
        self.gens = np.array([Gen(plane=plane) for i in range(CHROMOSOME_SIZE)])


class Genotype:
    """ Definition of Genotype - collection of Chromosomes """
    def __init__(self) -> None:
        self.chromosome_h = Chromosome(plane=Plane.HORIZONTAL)
        self.chromosome_v = Chromosome(plane=Plane.VERTICAL)

    def initialize(self, grid: Dict[str, np.array]) -> None:
        """
        Create random genotype which meets fallowing conditions:
            - Route start point and end point is based on defined START_POINT and END_POINT
            - Each arc/segment is tangent to previous and next arc/segment
        Params:                                                                     type:
        :param grid: holds x axis list, y axis list and mesh z values               Dict[str, np.array]
        :return: None
        """

        """ Initialize horizontal route path """
        self.__init_horizontal(grid=grid)
        route_length = self.__eval_h_route_len()

        """ Initialize vertical route path """
        self.__init_vertical(route_length=route_length)
        self.plot_init()

    def plot_init(self) -> None:
        if PLOT_INIT:
            plot_route_2d(plane=Plane.HORIZONTAL, p_dicts=[{'x': gen.point['x'], 'y': gen.point['y']} for gen in
                                                           self.chromosome_h.gens])
            plot_route_2d(plane=Plane.VERTICAL, p_dicts=[{'z': gen.point['z'], 'd': gen.point['d']} for gen in
                                                         self.chromosome_v.gens], init_tangent=0)

    def __init_horizontal(self, grid: Dict[str, np.array], sort_ascending=True) -> None:
        """
        Initialize Gens (Points) for hotizontal movement
        Params:                                                                     type:
        :param grid: holds x axis list, y axis list and mesh z values               Dict[str, np.array]
        :return: None
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
                x_drawn = (np.random.choice([i for i in range(MAP_LIMIT['xmin'] + DIST10KM,
                                                              MAP_LIMIT['xmax'] - DIST10KM)], 1))[0]
                y_drawn = (np.random.choice([i for i in range(MAP_LIMIT['ymin'] + DIST10KM,
                                                              MAP_LIMIT['ymax'] - DIST10KM)], 1))[0]

                """ Check if selected point has proper z value """
                if is_point_in_range(x_req=x_drawn, y_req=y_drawn, grid=grid):
                    gens[gen_id].point['x'] = x_drawn
                    gens[gen_id].point['y'] = y_drawn
                    break
                else:
                    """" Drawn point not in range - draw another one """
                    pass

        if sort_ascending:
            values_x = deepcopy([gen.point['x'] for gen in gens[1:-1]])
            values_y = deepcopy([gen.point['y'] for gen in gens[1:-1]])
            values_x.sort()
            values_y.sort()

            for gen_id in range(1, len(gens) - 1):
                (gens[gen_id]).point['x'] = values_x[gen_id - 1]
                (gens[gen_id]).point['y'] = values_y[gen_id - 1]

        print('Genotype: initialized horizontal: {}'.format([gen.point for gen in gens]))

    def __eval_h_route_len(self) -> float:
        """
        Use HypGeo module to evuate route length in horizontal
        Params:                                                                     type:
        :return: route length in horizontal plane                                   float
        """
        return evaluate_route_len(plane=Plane.HORIZONTAL, p_dicts=[{'x': gen.point['x'], 'y': gen.point['y']}
                                                                   for gen in self.chromosome_h.gens])

    def __init_vertical(self, route_length: float) -> None:
        """
        Initialize Gens (Points) for vertical movement
        Params:                                                                     type:
        :param route_length: route length in horizontal plane                       float
        :return: None
        """

        """ Initialize start point and end point coordinates """
        gens = self.chromosome_v.gens
        gens[0].point['z'] = START_POINT['z']
        gens[0].point['d'] = 0
        gens[-1].point['z'] = END_POINT['z']
        gens[-1].point['d'] = route_length

        """ Draw intermediate points coordinates """
        z_drawns = np.random.choice([i for i in range(ROUTE_MIN_HIGHT, ROUTE_MAX_HIGHT)], CHROMOSOME_SIZE - 2)
        d_drawns = np.random.choice([i for i in range(1, math.floor(route_length))], CHROMOSOME_SIZE - 2)

        """ Check if d values are unique """
        seen = set()
        if any(i in seen or seen.add(i) for i in d_drawns):
            raise ValueError('Values in d_drawns not unique: {}!'.format(d_drawns))

        """ Sort values by growing distance d """
        df_vert = pd.DataFrame({'z': z_drawns, 'd': d_drawns})
        df_vert = df_vert.sort_values(by=['d'], ignore_index=True)

        """ Set intermediate points """
        for gen_id in range(1, CHROMOSOME_SIZE - 1):
            gens[gen_id].point['z'] = df_vert['z'][gen_id - 1]
            gens[gen_id].point['d'] = df_vert['d'][gen_id - 1]

        print('Genotype: initialized vertical: {}'.format([gen.point for gen in gens]))

    def init_predefined(self, data_vector: np.array) -> None:
        self.chromosome_h.gens[0].point['x'] = START_POINT['x']
        self.chromosome_h.gens[0].point['y'] = START_POINT['y']
        self.chromosome_h.gens[-1].point['x'] = END_POINT['x']
        self.chromosome_h.gens[-1].point['y'] = END_POINT['y']

        for gen_id in range(1, CHROMOSOME_SIZE - 1):
            self.chromosome_h.gens[gen_id].point['x'] = data_vector[2 * (gen_id - 1)]
            self.chromosome_h.gens[gen_id].point['y'] = data_vector[2 * (gen_id - 1) + 1]

        horizontal_route_length = self.__eval_h_route_len()

        self.chromosome_v.gens[0].point['z'] = START_POINT['z']
        self.chromosome_v.gens[0].point['d'] = 0
        self.chromosome_v.gens[-1].point['z'] = END_POINT['z']
        self.chromosome_v.gens[-1].point['d'] = horizontal_route_length

        for gen_id in range(1, CHROMOSOME_SIZE - 1):
            self.chromosome_v.gens[gen_id].point['d'] = data_vector[2 * (gen_id - 1) + int(len(data_vector) / 2)]
            self.chromosome_v.gens[gen_id].point['z'] = data_vector[2 * (gen_id - 1) + 1 + int(len(data_vector) / 2)]

    def __eval_v_route_len(self) -> float:
        """
        Use HypGeo module to evuate route length in vertical
        Params:                                                                     type:
        :return: route length in vertical plane                                     float
        """
        return evaluate_route_len(plane=Plane.VERTICAL, p_dicts=[{'z': gen.point['z'], 'd': gen.point['d']}
                                                                 for gen in self.chromosome_v.gens])

    def get_route_length(self) -> float:
        """
        Get whole route length
        Params:                                                                     type:
        :return: route length                                                       float
        """
        return self.__eval_v_route_len()

    def get_route_length_horizontal(self) -> float:
        return self.__eval_h_route_len()

    def calculate_cost(self) -> float:
        """ temporary cost calculation """
        return self.get_route_length()


class Population:
    """ Definition of Population - collection of Genotypes """
    def __init__(self, pop_size=POPULATION_SIZE) -> None:
        self.genotypes = np.array([Genotype() for i in range(pop_size)])

    def initialize(self, grid: Dict[str, np.array]) -> None:
        """
        Params:                                                                     type:
        :param grid: holds x axis list, y axis list and mesh z values               Dict[str, np.array]
        :return: None
        """
        for genotype in self.genotypes:
            print('Population: initialize genotype {}'.format((list(self.genotypes)).index(genotype)))
            genotype.initialize(grid=grid)

    def create_new_generation(self) -> None:

        """ Create new offspring """
        offspring = np.array([Genotype() for i in range(len(self.genotypes))])

        """ select 30 best parents """
        best_parents = self.__select_best_individuals(num_of_ind=NEWPOP_BEST_PARENTS_NUM)
        offspring[NEWPOP_BEST_PARENTS_START:NEWPOP_BEST_PARENTS_NUM] = best_parents

        """ use 15 best parents to cross over offspring and generate 50 of them """
        children = self.__cross_over(num_of_ind=NEWPOP_CHILDREN_NUM, parents=best_parents)
        offspring[NEWPOP_CHILDREN_START:NEWPOP_CHILDREN_START + NEWPOP_CHILDREN_NUM] = children

        """ create random 20 individuals """
        randoms = np.array([Genotype() for i in range(NEWPOP_RANDOM_NUM)])
        x, y, z = read_filtered_data()
        grid = {'x': x, 'y': y, 'z': z}
        for random_ind in randoms:
            random_ind.initialize(grid=grid)

    def __select_best_individuals(self, num_of_ind: int) -> np.array:
        costs_desc = [{'index': index, 'cost': genotype.calculate_cost()} for
                      (index, genotype) in zip(range(len(self.genotypes)), self.genotypes)]

        sorted_indexes = [index for cost, index in sorted(zip([item['cost'] for item in costs_desc],
                                                              [item['index'] for item in costs_desc]))]

        self.genotypes = np.array([self.genotypes[index] for index in sorted_indexes])
        return np.array(self.genotypes[0: num_of_ind])

    def __cross_over(self, num_of_ind: int, parents: np.array) -> np.array:
        children = np.array([Genotype() for i in range(num_of_ind)])

        for child in children:
            cross_mask = np.random.choice([0, 1], CHROMOSOME_SIZE)
            for i in range(CHROMOSOME_SIZE):
                """ cross horizontal """
                parent1_id, parent2_id = CROSS_SELECT[i]
                parent1, parent2 = parents[parent1_id], parents[parent2_id]
                child.chromosome_h.gens[i] = parent1.chromosome_h.gens[i] if cross_mask[i] == 0\
                    else parent2.chromosome_h.gens[i]

            cross_mask = np.random.choice([0, 1], CHROMOSOME_SIZE)
            for i in range(CHROMOSOME_SIZE):
                """ cross vertical """
                parent1_id, parent2_id = CROSS_SELECT[i]
                parent1, parent2 = parents[parent1_id], parents[parent2_id]
                child.chromosome_v.gens[i] = parent1.chromosome_v.gens[i] if cross_mask[i] == 0 \
                    else parent2.chromosome_v.gens[i]

        return children




def is_point_in_range(x_req: int, y_req: int, grid: Dict[str, np.array]) -> bool:
    """
    Check if drawn coordinates are in map range and if corrensponding z value in not np.inf
    Params:                                                                     type:
    :param x_req: Drawn x coordinate                                            int
    :param y_req: Drawn y coordinate                                            int
    :param grid: Holds x axis list, y axis list and mesh z values               Dict[str, np.array]
    :return: True if drawn point is in map range, else False                    bool
    """
    x_axis = grid['x']
    y_axis = grid['y']

    if (int(x_req) not in range((grid['x'])[0], (grid['x'])[-1])) or\
            (int(y_req) not in range((grid['y'])[0], (grid['y'])[-1])):
        raise ValueError('Requested point ({},{}) out of map range! Ranges: x:({},{}) y:({},{})'.
                         format(x_req, y_req, (grid['x'])[0], (grid['x'])[-1], (grid['y'])[0], (grid['y'])[-1]))

    x_id_high, x_id_low = np.inf, np.inf
    y_id_high, y_id_low = np.inf, np.inf
    for x_id in range(len(x_axis)):
        if x_axis[x_id] == x_req:
            x_id_low = x_id
            x_id_high = x_id
            break
        elif x_axis[x_id] > x_req:
            x_id_low = x_id - 1
            x_id_high = x_id
            break
        else:
            """ Keep iterating """
            pass

    for y_id in range(len(y_axis)):
        if y_axis[y_id] == y_req:
            y_id_low = y_id
            y_id_high = y_id
            break
        elif y_axis[y_id] > y_req:
            y_id_low = y_id
            y_id_high = y_id + 1
            break
        else:
            """ Keep iterating """
            pass

    if np.inf in[x_id_low, x_id_high, y_id_low, y_id_high]:
        raise ValueError('Could not place requested point ({},{}) on map!'.format(x_req, y_req))

    return ((grid['z'])[x_id_low][y_id_low] != np.inf and
            (grid['z'])[x_id_high][y_id_low] != np.inf and
            (grid['z'])[x_id_high][y_id_high] != np.inf and
            (grid['z'])[x_id_low][y_id_high] != np.inf)
