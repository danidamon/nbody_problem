import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def update(iteration):
    points.set_data(dfs[iteration][0],dfs[iteration][1])
    plt.draw()
    return points

if len(sys.argv) < 2:
    sys.exit("Error: Need number of particles")
elif len(sys.argv) > 2:
    sys.exit("Error: Too many arguments given")

nParticles = int(sys.argv[1])

data = pd.read_csv('particles.txt', sep=" ", header=None)

# get the number of dataframes to split into 
# each dataframe contains the particles from a single iteration'''
rows = data.shape[0]
splits = rows/nParticles
dfs = np.array_split(data, splits)

# plot each dataframe
fig, ax = plt.subplots()
points, = ax.plot([], [], 'bo')
plt.xlim(0,1)
plt.ylim(0,1)

ani = animation.FuncAnimation(fig, update, int(splits), interval=300, repeat=True, blit=False)

plt.show()

