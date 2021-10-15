import multiprocessing
import sim_2pgs
from tqdm import trange
import numpy as np
from multiprocessing import shared_memory, Pool


class SharableListRepr(shared_memory.ShareableList):
    def __repr__(self):
        return f'{list(self)}'


def create_lattice(sparsity):
    freshLattice = np.random.randint(0, sparsity, (L, L))  # the lattice randomly selected

    minFunc = lambda x: min(x, 2)
    arrayMin = np.vectorize(minFunc)
    return arrayMin(freshLattice)


def gen_iter():
    for i in trange(len(e_range), leave=False):
        for j in trange(sample_size, leave=False):
            yield i, j


def evolution(division, iterations, lattice, selectionRate, reproductionRate, colonisationRate, index, baseRate, cost,
              good1_per_cell, good2_per_cell, e_range, sample_size, defe, part_coop, coop):
    for k in range(division):
        for l in range(iterations // division):
            sim_2pgs.simulation_one_generation(lattice, selectionRate, reproductionRate, colonisationRate,
                                               e_range[index], baseRate, cost, good1_per_cell, good2_per_cell)
        density = sim_2pgs.calc_density(lattice)
        if density[0] == 0:
            defe[index] += 1 / sample_size
            break
        if density[1] == 0 or density[2]==0:
            part_coop[index] += 1 / sample_size
            break
        if density[3] == 0:
            coop[index] += 1/sample_size
            break


if __name__ == '__main__':
    L = 100  # number of rows/columns
    N = L ** 2  # number of sites in total
    selectionRate = 1
    reproductionRate = 1
    colonisationRate = 0
    baseRate = 0.2
    good1_per_cell = 0.1
    good2_per_cell = 0.1
    cost = [0, 0.005, 0.005, 0.01, 0]
    division = 1
    iterations = 1000
    deathRate = 0
    e_values = 5
    sample_size = 10
    sparsity = 10

    e_range = np.logspace(-3, 1, num=e_values, endpoint=True)
    defector_list = SharableListRepr([0 for i in range(e_values)])
    cooperator_list = SharableListRepr([0 for i in range(e_values)])
    part_cooperator_list = SharableListRepr([0 for i in range(e_values)])

    with Pool(multiprocessing.cpu_count()-1) as p:
        for i, j in gen_iter():
            p.apply_async(evolution, (division, iterations, create_lattice(sparsity), selectionRate, reproductionRate,
                                      colonisationRate, i, baseRate, cost, good1_per_cell, good2_per_cell, e_range,
                                      sample_size, defector_list, part_cooperator_list, cooperator_list))
        p.close()
        p.join()
    print(defector_list, flush=True)
    print(part_cooperator_list, flush=True)
    print(cooperator_list, flush=True)
