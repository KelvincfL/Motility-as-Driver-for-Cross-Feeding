import sim_2pgs
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.animation as animation
from tqdm import trange

### Parameters: Initialising values for simulation



L = 100  # number of rows/columns
N = L ** 2  # number of sites in total
selectionRate = 1
reproductionRate = 1
exchangeRate = 0.1
colonisationRate = 0
baseRate = 0.05
deathRate = 0.05
good1_per_cell = [0,0.3,0,0]
good2_per_cell = [0,0,0.3,0]
cost = [0,0.05,0.05,0]
sparsity = 100

fileName = f"Base {baseRate} D {deathRate} C1 {cost[1]} C2 {cost[2]} G1 {good1_per_cell[1]} G2 {good2_per_cell[2]} Sp {sparsity} test 1"

print(fileName)

freshLattice = np.random.randint(0, sparsity, (L, L))  # the lattice randomly selected

minFunc = lambda x: min(x,3)
arrayMin = np.vectorize(minFunc)
lattice = arrayMin(freshLattice)

### Initialising for plot

fps = 30
seconds = 30
snapshots = [lattice.copy()]

### Simulating
print('Simulating:')

def count_cells_in_lattice(matrix):
    count = [0,0,0,0]
    for row in matrix:
        for cell in row:
            count[cell] += 1
    
    return count

frequency = []

frequency.append(count_cells_in_lattice(lattice))

for i in trange(fps*seconds,leave=False):
    sim_2pgs.simulation_one_generation(lattice, selectionRate, reproductionRate, colonisationRate, exchangeRate, baseRate, deathRate, cost, good1_per_cell, good2_per_cell)

    cell_count = count_cells_in_lattice(lattice)
    frequency.append(cell_count)    

    snapshots.append(lattice.copy())

    if cell_count.count(0) == 2:
        break #extinction of all other cells
    

print('')
print('Generations: {}'.format(i))
### Animation of process

print('')
print('Converting:')

colors = ['red', 'yellow', 'blue', 'grey']
cmap = mpl.colors.ListedColormap(colors)

fig = plt.figure( figsize=(5,5))

a = snapshots[0]

im = plt.imshow(a, interpolation='nearest', extent=[0.5, 0.5+L, 0.5, 0.5+L], cmap=cmap,aspect='auto')

plt.axis('off')

def animate_func(i):

    #Indication that the function is working
    if i % fps == 0:
        print( '.', end ='' )

    im.set_array(snapshots[i])

    return [im]

anim = animation.FuncAnimation(fig, animate_func, frames = len(frequency), interval = 1000 / fps)

anim.save('{}.html'.format(fileName), fps=fps, extra_args=['-vcodec', 'libx264'])

plt.clf()

### Plot graph
Y = np.array(frequency)
X = np.array(range(len(frequency)))

plt.plot(X,Y.T[0], label='(0,0)', color='red')
plt.plot(X,Y.T[1], label='(1,0)', color='orange')
plt.plot(X,Y.T[2], label='(0,1)', color='blue')
plt.title('Frequencies')
plt.legend()
plt.show()