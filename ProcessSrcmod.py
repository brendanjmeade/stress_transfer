import numpy
import scipy
import pyproj

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

# load file and extract geometric coordiantes and slip distribution
srcmod_eventname = 's1999HECTOR01SALI'
F = scipy.io.loadmat(srcmod_eventname + '.mat')
F = F[srcmod_eventname]

# extract the coordiantes for the rectangular fault patches
x = F[0]['seg1geoX']
y = F[0]['seg1geoY']
z = F[0]['seg1geoZ']


# plot the rectangular fault patches
fig = plt.figure()
ax = fig.gca(projection='3d')
theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)
z = np.linspace(-2, 2, 100)
r = z**2 + 1
#x = r * np.sin(theta)
#y = r * np.cos(theta)
ax.plot(x[0], y[0], z[0])

plt.show()


# extract the slip distribution

