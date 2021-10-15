import sim_ratio
import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange

### Parameters

L = 100  # number of rows/columns
N = L ** 2  # number of sites in total
selectionRate = 1
reproductionRate = 1
colonisationRate = 0
baseRate = 0.2
good_per_cooperator = 0.1
cost = [0, 0.05, 0]

division = 10
sample_size = 10
iterations = 1000

sparsity = 10
ratio_range = [0.1,0.5,1,2,10]
defector_list = [0 for i in range(len(ratio_range))]
cooperator_list = [0 for i in range(len(ratio_range))]
cooperatorExchangeRate = 1


def calc_density(lattice):
    L = len(lattice[0])
    density = [0, 0, 0]
    for i in range(L):
        for j in range(L):
            density[lattice[i][j]] += 1
    return density

def create_lattice(sparsity):
    freshLattice = np.random.randint(0, sparsity, (L, L))  # the lattice randomly selected

    minFunc = lambda x: min(x,2)
    arrayMin = np.vectorize(minFunc)
    return arrayMin(freshLattice)

def gen_iter():
    for i in trange(len(ratio_range), leave=False):
        for j in trange(sample_size, leave=False):
            lattice = create_lattice(sparsity)
            yield i, j, lattice


for i, j, lattice in gen_iter():
    defectorExchangeRate = cooperatorExchangeRate * ratio_range[i]
    for k in range(division):
        for l in range(iterations//division):
            sim_ratio.simulation_one_generation(lattice, selectionRate, reproductionRate, colonisationRate, \
                                                baseRate, cost, good_per_cooperator, \
                                                defectorExchangeRate, cooperatorExchangeRate)
        density = sim_ratio.calc_density(lattice)

        if density[0] == 0:
            defector_list[i] += 1/sample_size
            break
        if density[1] == 0:
            cooperator_list[i] += 1/sample_size
            break

print(defector_list, flush=True)
print(cooperator_list,flush=True)
# plt.plot(ratio_range, defector_list, marker='o', color='blue', label='Defector')
# plt.plot(ratio_range, cooperator_list, marker='o', color='orange', label='Cooperator')
# plt.xlabel('Exchange Rate')
# plt.ylabel('Extinction Probability')
# plt.legend()
# plt.show()