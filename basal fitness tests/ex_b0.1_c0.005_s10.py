import sim_fitness
import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange

### Parameters

L = 100  # number of rows/columns
N = L ** 2  # number of sites in total
selectionRate = 1
reproductionRate = 1
colonisationRate = 0
baseRate = 0.1
good_per_cooperator = 0.1
cost = [0,0.005,0]
sparsity = 10


#Run times
division = 20
sample_size = 10
max_generations = 1500
e_values = 11

e_range = np.logspace(-5,1,num=e_values,endpoint=True, base=10)
defector_list = [0 for i in range(len(e_range))]
cooperator_list = [0 for i in range(len(e_range))]

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
    for i in trange(len(e_range), leave=False):
        for j in trange(sample_size, leave=False):
            lattice = create_lattice(sparsity)
            yield i, j, lattice


for i, j, lattice in gen_iter():
    for k in range(division):
        for l in range(max_generations//division):
            sim_fitness.simulation_one_generation(lattice, selectionRate, reproductionRate, colonisationRate, e_range[i], baseRate, cost, good_per_cooperator)
        density = sim_fitness.calc_density(lattice)
        
        if density[0] == 0:
            defector_list[i] += 1/sample_size
            break
        if density[1] == 0:
            cooperator_list[i] += 1/sample_size
            break

print('Cost: {}; Sparsity: {}; Base: {}'.format(cost[1],sparsity, baseRate),flush=True)
print(e_range,flush=True)
print(defector_list,flush=True)
print(cooperator_list,flush=True)

'''
plt.plot(e_range, defector_list, marker='o', color='blue', label='Defector')
plt.plot(e_range, cooperator_list, marker='o', color='orange', label='Cooperator')
plt.xscale('log')
plt.xlabel('Exchange Rate')
plt.ylabel('Extinction Probability')
plt.legend()
plt.show()
'''
