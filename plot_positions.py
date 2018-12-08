# plot the energies
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np
import matplotlib.mlab as mlab
import math
# input file
filename = 'pos.dat'

# import data
data = np.loadtxt(filename)
cols = np.size(data,1)
# initial size of plot window
plt.figure(figsize=(8,6))

# plot
for i in range(1, cols-3):
    plt.plot(data[:,0], data[:,i],'-', label='Particle '+str(i), Linewidth=0.5)

plt.plot(data[:,0], data[:,cols-3],'b-', label='$\mu_X(t)$', Linewidth=1.5)
plt.plot(data[:,0], data[:,cols-2],'r-', label='$\mu_X(t) \pm \sigma_X(t)$', Linewidth=1.5)
plt.plot(data[:,0], data[:,cols-1],'r-')#, label='$\mu_X(t) + \sigma_X(t)$', Linewidth=1.5)
# labels
plt.xlabel('Time / [ms]', fontsize=20)
plt.ylabel('Positions / [nm]', fontsize=20)

# legend
plt.legend(loc='upper right')
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=12)

# axis limits
# tick fontsize
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.title('Position of '+str(cols-4)+ ' molecules, $\mu$ and $\mu \pm \sigma$')

# display the plot

plt.show()
