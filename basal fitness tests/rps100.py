import sim_original
import matplotlib.pyplot as plt
import numpy as np
from tqdm import trange

L = 100  # sites per row
N = L**2    # total number of sites
sample_size = 20
iterations = int(0.1*N)  # number of iterations why????
selectionRate = 1
reproductionRate = 1

M_range = np.logspace(-5,-3, num=9, endpoint=True, base=10)
probability_list = [0 for i in range(len(M_range))]
division = 20


def gen_iter():
    for i in trange(len(M_range), leave=False):
        exchangeRate = M_range[i]*N/2   # m=2*epsilon* N^(-1)
        for j in trange(sample_size, leave=False):
            lattice = np.random.randint(0, 4, (L, L))
            yield i, j, exchangeRate, lattice


for i, j, exchangeRate, lattice in gen_iter():
    for k in range(division):
        for l in range(iterations//division):
            sim_original.simulation_one_generation(lattice, selectionRate, reproductionRate, 0, exchangeRate)
        density = sim_original.calc_density(lattice)

        if density[0] == 0 or density[1] == 0 or density[2] == 0:
            probability_list[i] += 1/sample_size
            break


print(probability_list, flush=True)

'''
plt.plot(M_range, probability_list, marker='o')
plt.xscale('log')
plt.xlabel('Mobility')
plt.ylabel('Extinction Probability')
plt.show()
'''
