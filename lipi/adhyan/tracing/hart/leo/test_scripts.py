import os
import numpy as np
import matplotlib.pyplot as plt
import raytrace as rt
import numpy as np
import math
from raytrace import implane
import rtcore

freq_hz = 300.e6 # Hz, radio wave frequency
grid=(100,100)
rect=(-1, -1, 1, 1)
obs=(215, 0, 0)  # Position of the observer (the earth's coordinates) in SGI system
rsph = 10  # Radius of the integrarion sphere
niter = 1500  #The maximum iterations to perform using the algorithm

grid = (100,100)
nx = int(grid[0])
ny = int(grid[1])
# spacing between tracked rays
stepx = nx // 10   # 10
stepy = ny // 10   # 10


test1 = rt.implane(grid,
                   rect,
                   obs,
                   rsph,
                   freq=freq_hz,
                   mode='TbrIQUV',
                   trkparms=['pos'],
                   trknpmax=niter)
test1.package = '/home/rohit/sort/SORT/raytrace_leo_old_version/raytrace/py_raytr_threaded'

trkrays = []
for i in range(0, nx, stepx):
        for j in range(0, ny, stepy):
                    trkrays.append([i, j])
                   
test1.trace(niter)
traj = test1.pos


i=50
f,ax=plt.subplots(1,1)
ax.plot(traj[i,:,1],traj[i,:,2],'o-')
plt.show()

