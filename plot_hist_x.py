# plot the energies
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np

# input file
filename = 'x_hist.dat'

# import data
data = np.loadtxt(filename)
cols = np.size(data,1)
# initial size of plot window
plt.figure(figsize=(8,6))
scale = 5*1E-4

time = [10, 500, 3000, 20000]
h = []
# plot
for i in range(0, 4):
    plt.hist(data[:,i], bins=100, alpha=0.7, rwidth=0.85, label='Time [ms]:'+ str(time[i]*scale))
# labels
plt.xlabel('Position / [nm]', fontsize=20)
plt.ylabel(' #trajectories', fontsize=20)
#plt.text(-100, 600, r'$Going towards Gaussian dist.$')

# legend
plt.legend(loc='upper left')
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=12) 

# axis limits
# tick fontsize
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.title('Histogram of position of Brownian particle with respect to time')

# display the plot


plt.show()
