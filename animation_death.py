import sim_death
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.animation as animation

### Parameters: Initialising values for simulation

L = 100  # number of rows/columns
N = L ** 2  # number of sites in total
selectionRate = 1
reproductionRate = 0.8
exchangeRate = 1
colonisationRate = 0.2
deathRate = 0.15
baseRate = 0.2
good_per_cooperater = 0.1
cost = [0,0.05,0]
sparsity = 4

fileName = f'Ex {exchangeRate} Ba {baseRate} De {deathRate} Go {good_per_cooperater} Cost {cost[1]} Sp {sparsity} Col {colonisationRate} test 1'

freshLattice = np.random.randint(0, sparsity, (L, L))  # the lattice randomly selected

minFunc = lambda x: min(x,2)
arrayMin = np.vectorize(minFunc)
lattice = arrayMin(freshLattice)

### Initialising for plot

fps = 30
seconds = 30
snapshots = [lattice.copy()]

### Simulating
print(fileName)
print('Simulating:')

def count_cells_in_lattice(matrix):
    count = [0,0,0]
    for row in matrix:
        for cell in row:
            count[cell] += 1
    
    return count

frequency = []

frequency.append(count_cells_in_lattice(lattice))

i = 0
extinct = False

while i < seconds*fps - 1 and extinct == False:
    sim_death.simulation_one_generation(lattice, selectionRate, reproductionRate, colonisationRate, exchangeRate, baseRate, deathRate, cost, good_per_cooperater)

    cell_count = count_cells_in_lattice(lattice)
    frequency.append(cell_count)    

    if i % fps == 0:
        print('.',end ='')
    snapshots.append(lattice.copy())

    if cell_count.count(0) == 1:
        extinct = True
    
    i += 1

print('')
print('Generations: {}'.format(i))
### Animation of process

print('')
print('Converting:')

colors = ['blue', 'yellow', 'grey']
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

plt.plot(X,Y.T[0], label='Defector', color='blue')
plt.plot(X,Y.T[1], label='Cooperator', color='yellow')
plt.title('Frequencies')
plt.legend()
plt.show()
