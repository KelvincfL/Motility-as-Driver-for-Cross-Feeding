import random
import math
import numpy as np


### Types of cells
# Empty: 4
# (1,1): 3
# (0,0): 0
# (1,0): 1
# (0,1): 2

# periodic boundary
# each generation is where each individual has reacted on average once


### Fitness: Reproduction rate of cell proportional to fitness

def fitness_formula(cell_type: int, neighbours, baseRate, cost, good1_per_cell, good2_per_cell):
    fitness = baseRate - cost[cell_type] + \
              min(np.array(neighbours) @ np.array(good1_per_cell), np.array(neighbours) @ np.array(good2_per_cell))
    return max(0, fitness)


def fitness(index: tuple, lattice, baseRate, cost, good1_per_cell, good2_per_cell):
    L = len(lattice[0])
    neighbours = [0, 0, 0, 0]
    for i in range(-1, 2):
        for j in range(-1, 2):
            neighbours[lattice[(index[0] + i) % L][(index[1] + j) % L]] += 1

    return fitness_formula(lattice[index[0]][index[1]], neighbours, baseRate, cost, good1_per_cell, good2_per_cell)


### Reaction: Selection or Reproduction
def reaction_single(lattice, selectionRate, reproductionRate, colonisationRate, baseRate, cost, good1_per_cell,
                    good2_per_cell):
    L = len(lattice[0])

    # Choosing random pair of neighbours
    site_to_update = random.choices(range(0, L), k=2)
    position = random.choice(
        [[-1, 0], [1, 0], [0, -1], [0, 1]])  # position of neighbour with respect to the site, random selection
    neighbour_index = ((site_to_update[0] + position[0]) % L, (site_to_update[1] + position[1]) % L)
    site_neighbour = lattice[neighbour_index[0]][neighbour_index[1]]

    # Determining type of reaction
    reaction_type = random.random()  # random number between 0 and 1

    # Reproduction if we have compatible pair (i.e. first cell empty)
    if reaction_type < reproductionRate / (selectionRate + reproductionRate + colonisationRate):

        # Valid reproduction set
        if lattice[site_to_update[0]][site_to_update[1]] == 4 and site_neighbour != 4:

            # Reproduction dependent on fitness
            if random.random() < fitness(neighbour_index, lattice, baseRate, cost, good1_per_cell, good2_per_cell):
                lattice[site_to_update[0]][site_to_update[1]] = site_neighbour

    # Selection if we have compatible pair (i.e. first cell has lower fitness) - can kill or colonise
    elif lattice[site_to_update[0]][site_to_update[1]] != site_neighbour \
            and (site_neighbour != 4 and lattice[site_to_update[0]][site_to_update[1]] != 4):  # for selection

        relative_fitness = fitness(neighbour_index, lattice, baseRate, cost, good1_per_cell, good2_per_cell) \
                           - fitness(site_to_update, lattice, baseRate, cost, good1_per_cell,
                                     good2_per_cell)  # of neighbour vs update cell

        if relative_fitness > 0 or (
                relative_fitness == 0 and random.random() < 0.5):  # greater fitness wins or equal prob
            if reaction_type < (reproductionRate + selectionRate) / (
                    selectionRate + reproductionRate + colonisationRate):
                lattice[site_to_update[0]][site_to_update[1]] = 4  # kill
            else:
                lattice[site_to_update[0]][site_to_update[1]] = site_neighbour  # else colonise


### Diffusion by n_diff exchanges of nearest neighbours
def diffusion(lattice, n_diff):
    L = len(lattice[0])
    i = 0
    while i < n_diff:

        # Choose random pair of neighbours
        site_to_update = random.choices(range(0, L), k=2)
        position = random.choice(
            [[-1, 0], [1, 0], [0, -1], [0, 1]])  # position of neighbour with respect to the site, random selection
        site_neighbour = lattice[(site_to_update[0] + position[0]) % L][(site_to_update[1] + position[1]) % L]

        # Exchange (trivial if both cells empty)
        if site_neighbour != lattice[site_to_update[0]][site_to_update[1]]:
            lattice[(site_to_update[0] + position[0]) % L][(site_to_update[1] + position[1]) % L] = \
                lattice[site_to_update[0]][site_to_update[1]]
            lattice[site_to_update[0]][site_to_update[1]] = site_neighbour
        i += 1


'''For large epsilon (epsilon is the exchange rate), many diffusion events will take place 
between single reactions. Calculate the distribution of having n_diff diffusions between
two reaction events, =epsilon^n_diff*(1-epsilon)  (poisson distribution with unit time). 
The probability distribution of having n_diff<=m is
then 1-epsilon^(m+1). Draw n_diff randomly according to this
distribution: n_diff+1=log(r)/log(epsilon/(1+epsilon)) with r uniformly distributed random
number in [0,1]'''


def simulation_one_generation(lattice, selectionRate, reproductionRate, colonisationRate, exchangeRate, baseRate, cost,
                              good1_per_cell, good2_per_cell):
    L = len(lattice[0])
    reaction_rate = selectionRate + reproductionRate + exchangeRate + colonisationRate
    i = 0

    # One generation: each cell has 1 reaction on average
    while i < (L ** 2) * reaction_rate:
        # Geometric random variable: number of exchanges between reaction processes
        n_diff = abs(math.log(random.random()) / math.log(exchangeRate / reaction_rate))

        diffusion(lattice, n_diff)
        reaction_single(lattice, selectionRate, reproductionRate, colonisationRate, baseRate, cost, good1_per_cell,
                        good2_per_cell)
        i += 1


def calc_density(lattice):
    L = len(lattice[0])
    density = [0, 0, 0, 0, 0]
    for i in range(L):
        for j in range(L):
            density[lattice[i][j]] += 1
    return density
