#!/usr/bin/env python
"""
This script setup initial condition and parameters for 2D diffusion, call 2D diffusion numerical solver, and display solution
"""
import scipy as sp
import matplotlib.pyplot as plt
from pylab import *
import diffuse2Dsolver

# Setup parameters:
a=.1          # Diffusion constant.
timeTotal=.1  # Number of time-steps to evolve system.

# Set up initial conditions
dx=0.01        # Interval size in x-direction.
dy=0.01        # Interval size in y-direction.

# Now, set the initial conditions (ui). Here we just start from a filled
# circular shape of radius 0.001
nx = int(1/dx)
ny = int(1/dy)
ui = sp.zeros([nx,ny])
pp = ui+ui
for i in range(nx):
	for j in range(ny):
		if ( ( (i*dx-0.5)**2+(j*dy-0.5)**2 <= 0.001)):
				ui[i,j] = 0.1

# run solver of diffusion equation                                
outputMat = diffuse2Dsolver.diffSolve(ui, a, dx, dy, timeTotal)


#display result. We have 6 subfigure showing results from different diffusion time points. 
intv_ind = int(len(outputMat)/5) # interval index between each plot
intv_time = timeTotal/5 # interval time between each plot

f, axarr = plt.subplots(2, 3)
u = outputMat[0]
axarr[0, 0].imshow( u, cmap=cm.hot, interpolation='nearest', origin='lower')
axarr[0, 0].set_title('time = 0')
u = outputMat[intv_ind - 1]
axarr[0, 1].imshow(u, cmap=cm.hot, interpolation='nearest', origin='lower')
axarr[0, 1].set_title('time = ' + str(intv_time))
u = outputMat[2*intv_ind - 1]
axarr[0, 2].imshow(u, cmap=cm.hot, interpolation='nearest', origin='lower')
axarr[0, 2].set_title('time = ' + str(2*intv_time))
u = outputMat[3*intv_ind - 1]
axarr[1, 0].imshow(u, cmap=cm.hot, interpolation='nearest', origin='lower')
axarr[1, 0].set_title('time = ' + str(3*intv_time))
u = outputMat[4*intv_ind - 1]
axarr[1, 1].imshow(u, cmap=cm.hot, interpolation='nearest', origin='lower')
axarr[1, 1].set_title('time = ' + str(4*intv_time))
u = outputMat[5*intv_ind - 1]
axarr[1, 2].imshow(u, cmap=cm.hot, interpolation='nearest', origin='lower')
axarr[1, 2].set_title('time = ' + str(5*intv_time))

plt.show()
