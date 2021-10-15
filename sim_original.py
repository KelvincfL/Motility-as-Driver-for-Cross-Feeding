import random
import math

# periodic boundary
# each generation is where each individual has reacted on average once



# calculation
# reaction step between individuals on two neighbouring site, involves birth and death process

### Reaction: Selection or Reproduction
def reaction_single(lattice, selectionRate, reproductionRate, colonisationRate):
    L = len(lattice[0])

    #Choosing random pair of neighbours
    site_to_update = random.choices(range(0, L), k=2)
    position = random.choice(
        [[-1, 0], [1, 0], [0, -1], [0, 1]])  # position of neighbour with respect to the site, random selection
    site_neighbour = lattice[(site_to_update[0] + position[0]) % L][(site_to_update[1] + position[1]) % L]

    #Determining type of reaction
    reaction_type = random.random()  # random number between 0 and 1

    #Reproduction if we have compatible pair (i.e. first cell empty)
    if reaction_type < reproductionRate / (selectionRate + reproductionRate + colonisationRate):
        if lattice[site_to_update[0]][site_to_update[1]] == 3:
            lattice[site_to_update[0]][site_to_update[1]] = site_neighbour

    #Selection if we have compatible pair (i.e. first cell is weaker) - can kill or colonise
    elif lattice[site_to_update[0]][site_to_update[1]] == ((site_neighbour + 1) % 3) \
            and (site_neighbour != 3 and lattice[site_to_update[0]][site_to_update[1]] != 3):  # for selection
        if reaction_type < (reproductionRate + selectionRate) / (selectionRate + reproductionRate + colonisationRate):
            lattice[site_to_update[0]][site_to_update[1]] = 3 #kill
        else:
            lattice[site_to_update[0]][site_to_update[1]] = site_neighbour  # else colonise


### Diffusion by n_diff exchanges of nearest neighbours
def diffusion(lattice, n_diff):
    L = len(lattice[0])
    i = 0
    while i < n_diff:

        #Choose random pair of neighbours
        site_to_update = random.choices(range(0, L), k=2)
        position = random.choice(
            [[-1, 0], [1, 0], [0, -1], [0, 1]])  # position of neighbour with respect to the site, random selection
        site_neighbour = lattice[(site_to_update[0] + position[0]) % L][(site_to_update[1] + position[1]) % L]

        #Exchange (trivia if both cells empty)
        if site_neighbour != lattice[site_to_update[0]][site_to_update[1]]:
            lattice[(site_to_update[0] + position[0]) % L][(site_to_update[1] + position[1]) % L] = \
                lattice[site_to_update[0]][site_to_update[1]]
            lattice[site_to_update[0]][site_to_update[1]] = site_neighbour
        i += 1


def simulation_one_generation(lattice, selectionRate, reproductionRate, colonisationRate, exchangeRate):
    L = len(lattice[0])
    reaction_rate = selectionRate + reproductionRate + exchangeRate + colonisationRate
    i = 0

    #One generation: each cell has 1 reaction on average
    while i < (L ** 2) * reaction_rate:

        #Geometric random variable: number of exchanges between reaction processes
        n_diff = abs(math.log(random.random()) / math.log(exchangeRate / reaction_rate))

        diffusion(lattice, n_diff)
        reaction_single(lattice, selectionRate, reproductionRate, colonisationRate)
        i += 1

def calc_density(lattice):
    L = len(lattice[0])
    density = [0, 0, 0, 0]
    for i in range(L):
        for j in range(L):
            density[lattice[i][j]] += 1
    return density
