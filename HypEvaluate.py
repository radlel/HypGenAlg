from HypData import read_filtered_data
from HypGenAlg import Population
from HypTypes import *


# def fitness(data_vector: np.array) -> float:
#     """
#     Expected input in data_vector:
#     data_vector[0] =Point[1].x
#     data_vector[1] =Point[1].y
#     data_vector[2] =Point[2].x
#     data_vector[3] =Point[2].y
#     data_vector[4] =Point[3].x
#     data_vector[5] =Point[3].y
#     data_vector[6] =Point[1].d
#     data_vector[7] =Point[1].z
#     data_vector[8] =Point[2].d
#     data_vector[9] =Point[2].z
#     data_vector[10] =Point[3].d
#     data_vector[11] =Point[3].z
#     """
#     penalty = 0
#
#     """ Get grid in order to check if value is valid, if not penalty must be added """
#     x, y, z = read_filtered_data()
#     grid = {'x': x, 'y': y, 'z': z}
#
#     """ Check if point is in allowed range """
#     for index in [i for i in range(int(len(data_vector) / 2)) if i % 2 == 0]:
#         if not is_point_in_range(x_req=data_vector[index], y_req=data_vector[index + 1], grid=grid):
#             penalty += 10000000
#             print('penalty: {}'.format(penalty))
#
#     """ Initialize Genotype """
#     genotype = Genotype()
#     genotype.init_predefined(data_vector=data_vector)
#     route_len = genotype.get_route_length()
#     print('RL:', round(genotype.get_route_length()))
#     return route_len + penalty


if '__main__' == __name__:
    pop = Population(pop_size=1)
    x, y, z = read_filtered_data()
    grid = {'x': x, 'y': y, 'z': z}
    pop.initialize(grid=grid)
