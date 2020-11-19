from HypData import read_filtered_data
from HypGenAlg import Population


if '__main__' == __name__:
    pop = Population(pop_size=1)
    x, y, z = read_filtered_data()
    grid = {'x': x, 'y': y, 'z': z}
    pop.initialize(grid=grid)
