from HypGenAlg import Population, plot_route_2d

from HypTypes import *


if '__main__' == __name__:
    pop = Population()

    pop.initialize_random()
    for i in range(GENERATIONS_NUM):
        print('Create new generation', i)
        pop.create_new_generation()
        best_ind = pop.get_best_individual()
        print('>>> BEST IND FITNESS:', best_ind.fenotype.fitness())

    """ Plot best route """
    best_route_desc = best_ind.fenotype.get_route_desc()
    best_route_len = best_ind.fenotype.route_horizontal_length
    best_p_descs = best_ind.genotype.get_points_descs(plane=Plane.HORIZONTAL)
    plot_route_2d(plane=Plane.HORIZONTAL, route_desc=best_route_desc, route_len=best_route_len, p_dicts=best_p_descs)
