import sim_rps
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.animation as animation

###Initialising values for simulation

fileName = '0.00000000000001 - test 3'

L = 50  # number of rows/columns
N = L ** 2  # number of sites in total
selectionRate = 1
reproductionRate = 1
exchangeRate = 0.00000000000001
colonisationRate = 0
lattice = np.random.randint(0, 4, (L, L))  # the lattice randomly selected

###Initialising for plot

fps = 50
seconds = 50
snapshots = [lattice.copy()]

### Simulating matrix
print('Simulating:')

def count_cells_in_lattice(matrix):
    count = [0,0,0,0]
    for row in matrix:
        for cell in row:
            count[cell] += 1
    
    return count

frequency = []

frequency.append(count_cells_in_lattice(lattice))

i = 0
extinct = False

while i < seconds*fps - 1 and extinct == False:
    ### Simulate generation
    sim_rps.simulation_one_generation(lattice, selectionRate, reproductionRate, colonisationRate, exchangeRate)
    cell_count = count_cells_in_lattice(lattice)
    frequency.append(cell_count)

    if i % fps == 0:
        print('.',end ='')
    snapshots.append(lattice.copy())
    
    if cell_count.count(0) == 3:
        extinct = True
    i += 1

### Animation of process
print('')
print('Converting:')

colors = ['blue', 'red', 'yellow', 'grey']
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

plt.plot(X,Y.T[0], label='Rock', color='blue')
plt.plot(X,Y.T[1], label='Scissors', color='red')
plt.plot(X,Y.T[2], label='Paper', color='yellow')
plt.title('Frequencies')
plt.legend()
plt.show()
