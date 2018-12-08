# plot the energies
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np
import math
import matplotlib.mlab as mlab
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
timestep = 1E-7
radius = (2.79*1E-6)/2.0;
density = 2.65 * 1E3;
m = 4.0*math.pi*radius*radius*radius *density*3.0; # assuming particle is sphere-ish
temperature = 297;
tau = 48.5*1E-6; # 48.5 | 147.3;
eta = 1.0/tau;
c_0 = math.exp(-2*eta*timestep);
k_b = 1.380*1E-23;
v0 = 2
mu = v0*math.exp(-eta*timestep)

variance = k_b*temperature *(1-math.exp(-2*eta*timestep)) / m
sigma = math.sqrt(variance)*5*1E5
x = np.linspace(mu - 3*sigma, mu + 3*sigma, 1000)
plt.plot(x,mlab.normpdf(x, mu, sigma), linewidth = 2)

print(mu, sigma)


plt.show()
