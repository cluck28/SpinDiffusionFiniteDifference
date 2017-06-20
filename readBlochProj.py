import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
import matplotlib.cm as cm
import numpy as np
import math
from scipy import optimize

list = np.array([90])

item = float(list[0])

#Read in data
data = np.loadtxt('MagnitudeMagnetizationEcho_%dPulse.txt' % (item))

#Unpack data
gridSize = len(data[0,:])-1
times = len(data[:,0])
time = np.zeros((times,gridSize))
zpos = np.zeros((times,gridSize))
for i in range(times):
    for j in range(gridSize):
        time[i][j] = data[i,0]
for i in range(times):
    for j in range(gridSize):
        zpos[i][j] = j*1.0
data = data[:,1:gridSize+1]
fig = plt.figure()
ax = Axes3D(fig)
surf = ax.plot_wireframe(zpos,time,data,rstride=1000,cstride=10)


ax.set_xlabel('X')
#ax.set_xlim3d(0, 1000)
ax.set_ylabel('Y')
#ax.set_ylim3d(0, 0.005)
ax.set_zlabel('Z')
#ax.set_zlim3d(-1, 1)

plt.show()
