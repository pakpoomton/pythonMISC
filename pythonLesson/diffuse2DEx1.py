#!/usr/bin/env python
"""
A program which uses an explicit finite difference
scheme to solve the diffusion equation with fixed
boundary values and a given initial value for the
density.

Two steps of the solution are stored: the current
solution, u, and the previous step, ui. At each time-
step, u is calculated from ui. u is moved to ui at the
end of each time-step to move forward in time.

"""
import scipy as sp
import time
import matplotlib.pyplot as plt
from pylab import *
import diffuse2Dsolver

# Declare some variables:
dx=0.01        # Interval size in x-direction.
dy=0.01        # Interval size in y-direction.
a=.5          # Diffusion constant.



timesteps=250  # Number of time-steps to evolve system.


# Now, set the initial conditions (ui).
nx = int(1/dx)
ny = int(1/dy)
ui = sp.zeros([nx,ny])
pp = ui+ui
for i in range(nx):
	for j in range(ny):
		if ( ( (i*dx-0.5)**2+(j*dy-0.5)**2 <= 0.08)
			& ((i*dx-0.5)**2+(j*dy-0.5)**2>=.07) ):
				ui[i,j] = 1

# run solver of diffusion equation                                
outputMat = diffuse2Dsolver.diffSolve(ui, a, dx, dy, timesteps)

f1 = plt.figure()
u = outputMat[0]
# display diffusion output   
plt.imshow( u, cmap=cm.hot, interpolation='nearest', origin='lower')

f2 = plt.figure()
u = outputMat[5]
# display diffusion output   
plt.imshow( u, cmap=cm.hot, interpolation='nearest', origin='lower')


f3 = plt.figure()
u = outputMat[25]
# display diffusion output   
plt.imshow( u, cmap=cm.hot, interpolation='nearest', origin='lower')

f4 = plt.figure()
u = outputMat[125]
# display diffusion output   
plt.imshow( u, cmap=cm.hot, interpolation='nearest', origin='lower')


plt.show()
