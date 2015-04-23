import scipy as sp
import time
import matplotlib.pyplot as plt
from pylab import *


def evolve_ts(u, ui, a, dx, dy, dx2, dy2, dt):
   """
   This function uses a numpy expression to
   evaluate the derivatives in the Laplacian, and
   calculates u[i,j] based on ui[i,j].
   """
   u[1:-1, 1:-1] = ui[1:-1, 1:-1] + a*dt*( (ui[2:, 1:-1] - 2*ui[1:-1, 1:-1] + ui[:-2, 1:-1])/dx2 + (ui[1:-1, 2:] - 2*ui[1:-1, 1:-1] + ui[1:-1, :-2])/dy2 )


def diffSolve(ui, a, dx, dy , timesteps):
    dx2=dx**2 # To save CPU cycles, we'll compute Delta x^2
    dy2=dy**2 # and Delta y^2 only once and store them.

    # For stability, this is the largest interval possible
    # for the size of the time-step:
    dt = dx2*dy2/( 2*a*(dx2+dy2) )
    

    # Start u and ui off as zero matrices:
    nx = int(1/dx)
    ny = int(1/dy)
    u = sp.zeros([nx,ny])

    outputMat = []
    outputMat.append(ui)
    
    # Now, start the time evolution calculation...
    tstart = time.time()
    for m in range(1, timesteps+1):
	evolve_ts(u, ui, a, dx, dy, dx2, dy2, dt)
        ui = sp.copy(u)
        outputMat.append(sp.copy(u))
	print "Computing u for m =", m
    tfinish = time.time()
    print "Done."
    print "Total time: ", tfinish-tstart, "s"
    print "Average time per time-step using numpy: ", ( tfinish - tstart )/timesteps, "s."
    return outputMat
