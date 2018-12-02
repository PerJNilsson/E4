# plot the energies
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np
import math

# input file
filename = 'v_hist.dat'


radius = (2.79*1E-6)/2.0
density = 2.65 * 1E3
m = 4.0*math.pi*radius*radius*radius *density*3.0
k_b = 1.380*1E-23
T = 297


# import data
data = np.loadtxt(filename)
cols = np.size(data,1)
# initial size of plot window
plt.figure(figsize=(8,6))
scale = 5*1*1E-5
time = [10*scale, 500*scale, 3000*scale, 20000*scale]
h = []
# plot
for i in range(0, 4):
    plt.hist(data[:,i], bins=100, alpha=0.7, rwidth=0.85, label='Time [ms]:'+ str(round(time[i],5)))

# labels
plt.xlabel('Velocities / [mm/s]', fontsize=20)
plt.ylabel(' #trajectories', fontsize=20)

# legend
plt.legend(loc='upper left')
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=12) 

# axis limits
# tick fontsize
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.title('Historgram of velocities of a Brownian particle with repect to time')

# display the plot


plt.show()
