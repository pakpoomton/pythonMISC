import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.animation as animation

## This script shows how to solve 1D heat equation numerically.
# Here we shave start from a hot spot in the middle of a 1D rod

## initialise parameters of heat equation
depth = 1. # length of the rod
timeTotal = .2 # simulation time 
K = 1 # diffusion coefficient


## specify parameter for numerical solution
Nz = 40. # number of discrete step in space
dz = depth/Nz # size of discrete step in space
# array of each location in our discrete space
# size of discrete step in time, note that
# dt<=dz**2/2/K for stability of numerical solution
global dt
dt = 0.9*dz**2/4/K 
Z_range = np.arange(0,depth+dz,dz)
t_range = np.arange(0, timeTotal, dt)
Nt = int(timeTotal/dt) # total number of discrete time points
midpoint = int(Nz/2) # mid point in discrete space


## initialise 3D solution matrix
# each row&column for each different locations; each layer for each time point
global T
T = np.zeros((Nt+1,Nz+1, Nz+1)) 
T[0, midpoint, midpoint] = 100 # initial temperature

## iterate through each time point and calculate 
for i in range(1, Nt):

    #second derivative of temperature in space   
    Ttop = T[i-1,0:-2,1:-1]
    Tleft = T[i-1,1:-1,0:-2]
    Tbottom = T[i-1,2:,1:-1]
    Tright = T[i-1,1:-1,2:]
    Tcenter = T[i-1,1:-1,1:-1]
    depth_2D = (Ttop+Tleft+Tbottom+Tright-4*Tcenter)/dz**2
    #first derivative of temperature in time
    time_1D = K*depth_2D
    #calculate next time point temperature for each non-boundary point
    T[i,1:-1,1:-1] = Tcenter+dt*time_1D
    #apply boundary condition on both ends saying that no heat going in/out
    T[i,0,:] = T[i,1,:]
    T[i,-1,:] = T[i,-2,:]
    T[i,:,0] = T[i,:,1]
    T[i,:,-1] = T[i,:,-2]



## display the result
zStep = int(Nz/4) # plot five locations on x-axis
xbound = np.arange(0, Nz+1, zStep)
xlabel = map(str, xbound*depth/Nz)

sampleT= int(Nt/24)
Tshow = T[sampleT:,:,:]
Tmin = Tshow.min()
Tmax = Tshow.max()

fig, ax = plt.subplots()


def animate(t):
    im = plt.imshow(T[t,:,:], cmap=plt.cm.hot, vmin=Tmin, vmax=Tmax)
    plt.title('time = '+str(round(float(t*dt),3)))
    plt.xticks(xbound, xlabel)
    plt.yticks(xbound, xlabel) 

 
ani = animation.FuncAnimation(fig, animate, np.arange(1,Nt,4),interval=1)

plt.colorbar(orientation='vertical',fraction=0.046, pad=0.04)
plt.show()








'''
# we will plot heat distribution from six different time points
sampleT= int(Nt/6)
sampleTarray = sampleT*np.array([1,2,3,4,5,6])*timeTotal/Nt
sampleTarray = np.round(sampleTarray,2)
sampleTarray = map(str,sampleTarray)

Tshow =np.array([T[1*sampleT,:,:], T[2*sampleT,:,:], T[3*sampleT,:,:], T[4*sampleT,:,:], T[5*sampleT,:,:], T[6*sampleT,:,:]])

Tmin = Tshow.min()
Tmax = Tshow.max()


zStep = int(Nz/4) # plot five locations on x-axis
xbound = np.arange(0, Nz+1, zStep)
xlabel = map(str, xbound*depth/Nz)


fig, axes = plt.subplots(nrows=2, ncols=3)
i = 0
for ax in axes.flat:
    im = ax.imshow(Tshow[i], cmap=plt.cm.hot, vmin=Tmin, vmax=Tmax)
    ax.set_xticks(xbound)
    ax.set_xticklabels(xlabel) 
    ax.set_yticks(xbound)
    ax.set_yticklabels(xlabel) 
    ax.set_title('t = ' + sampleTarray[i])
    i = i+1

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax=cbar_ax)

plt.show()
'''

