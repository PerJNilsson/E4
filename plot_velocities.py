# plot the energies
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np

# input file
filename = 'velocities.dat'

# import data
data = np.loadtxt(filename)
cols = np.size(data,1)

# initial size of plot window
plt.figure(figsize=(8,6))

# plot
for i in range(1, cols-3):
    plt.plot(data[:,0], data[:,i],'-', label='Trajectory '+str(i), linewidth=0.5)
plt.plot(data[:,0], data[:,cols-3],'b-', label='$\mu_v(t)$')
plt.plot(data[:,0], data[:,cols-2],'r-', label='$\mu_v(t) \pm \sigma_v(t)$', Linewidth=1.5)
plt.plot(data[:,0], data[:,cols-1],'r-')#, label='$\mu + \sigma$')

# labels
plt.xlabel('Time / [ms]', fontsize=20)
plt.ylabel('Velocity / [mm/s]', fontsize=20)

# legend
plt.legend(loc='upper right')
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=12) 

# axis limits
# tick fontsize
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.title('Velocities of particles')

# display the plot


plt.show()
