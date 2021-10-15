import random
import math
import numpy as np

### Types of cells
# Empty: 2
# (0): 0
# (1): 1

# periodic boundary
# each generation is where each individual has reacted on average once


### Fitness: Reproduction rate of cell proportional to fitness

def fitness_formula(cell_type: int, cooperators, defectors, baseRate, cost, good_per_cooperator):
    fitness = baseRate - cost[cell_type] + cooperators * good_per_cooperator
    return max(0, fitness)


def fitness(index: tuple, lattice, baseRate, cost, good_per_cooperator):
    L = len(lattice[0])
    neighbours = [0, 0, 0]
    for i in range(-1, 2):
        for j in range(-1, 2):
            neighbours[lattice[(index[0] + i) % L][(index[1] + j) % L]] += 1

    return fitness_formula(lattice[index[0]][index[1]], neighbours[1], neighbours[0], baseRate, cost,
                           good_per_cooperator)


### Reaction: Selection or Reproduction
def reaction_single(lattice, selectionRate, reproductionRate, colonisationRate, baseRate, cost, good_per_cooperator):
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
        if lattice[site_to_update[0]][site_to_update[1]] == 2 and site_neighbour != 2:

            # Reproduction dependent on fitness
            if random.random() < fitness(neighbour_index, lattice, baseRate, cost, good_per_cooperator):
                lattice[site_to_update[0]][site_to_update[1]] = site_neighbour

    # Selection if we have compatible pair (i.e. first cell has lower fitness) - can kill or colonise
    elif lattice[site_to_update[0]][site_to_update[1]] != site_neighbour \
            and (site_neighbour != 2 and lattice[site_to_update[0]][site_to_update[1]] != 2):  # for selection

        relative_fitness = fitness(neighbour_index, lattice, baseRate, cost, good_per_cooperator) \
                           - fitness(site_to_update, lattice, baseRate, cost,
                                     good_per_cooperator)  # of neighbour vs update cell

        if relative_fitness > 0 or (
                relative_fitness == 0 and random.random() < 0.5):  # greater fitness wins or equal prob
            if reaction_type < (reproductionRate + selectionRate) / (
                    selectionRate + reproductionRate + colonisationRate):
                lattice[site_to_update[0]][site_to_update[1]] = 2  # kill
            else:
                lattice[site_to_update[0]][site_to_update[1]] = site_neighbour  # else colonise


# Diffusion by n_diff exchanges of nearest neighbours
# ratio should be interpreted as defector/cooperator
def diffusion(lattice, n_diff, baseRate, cost, good_per_cooperator, bias_factor):
    L = len(lattice[0])
    i = 0
    while i < n_diff:
        # Choose random site
        site_to_update = random.choices(range(0, L), k=2)
        possibles = [[-1, 0], [1, 0], [0, -1], [0, 1]]
        fitnesses = [0, 0, 0, 0]
        for i in range(len(possibles)):
            fitnesses[i] = (fitness(((site_to_update[0] + possibles[i][0]) % L,\
                            (site_to_update[1] + possibles[i][1]) % L), lattice, baseRate, cost, good_per_cooperator))**bias_factor
        temp_total = sum(fitnesses)
        # at least one neighbour has fitness >0
        if temp_total > 0:
            # position of neighbour with respect to the site, weighted random selection
            position = random.choices(possibles, weights=fitnesses)
        else:
            position = random.choice(possibles)
        site_neighbour = lattice[(site_to_update[0] + position[0][0]) % L][(site_to_update[1] + position[0][1]) % L]

        # Exchange (trivial if both cells empty)
        if site_neighbour != lattice[site_to_update[0]][site_to_update[1]]:
            lattice[(site_to_update[0] + position[0][0]) % L][(site_to_update[1] + position[0][1]) % L] = \
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


def simulation_one_generation(lattice, selectionRate, reproductionRate, colonisationRate, baseRate, cost,
                              good_per_cooperator, exchangeRate, bias_factor):
    L = len(lattice[0])
    reaction_rate = selectionRate + reproductionRate + exchangeRate + colonisationRate
    i = 0

    # One generation: each cell has 1 reaction on average
    while i < (L ** 2) * reaction_rate:
        # Geometric random variable: number of exchanges between reaction processes
        n_diff = abs(math.log(random.random()) / math.log(exchangeRate / reaction_rate))

        diffusion(lattice, n_diff, baseRate, cost, good_per_cooperator, bias_factor)
        reaction_single(lattice, selectionRate, reproductionRate, colonisationRate, baseRate, cost, good_per_cooperator)
        i += 1


def calc_density(lattice):
    L = len(lattice[0])
    density = [0, 0, 0]
    for i in range(L):
        for j in range(L):
            density[lattice[i][j]] += 1
    return density
