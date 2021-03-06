import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.animation as animation

## This script shows how to solve 2D reaction-diffusion numerically.

## initialise parameters of heat equation
depth = 100 # length of the rod
timeTotal = 4.0 # simulation time 



k_u = 1 # diffusion rate of u
k_v = 40 # diffusion rate of v
c1 = 0.1
c2 = 0.9
cm1 = 1
c3 = 1


## specify parameter for numerical solution
Nz = 400. # number of discrete step in space
dz = depth/Nz # size of discrete step in space
# array of each location in our discrete space
# size of discrete step in time, note that
# dt<=dz**2/4/K for stability of numerical solution
global dt
dt = 0.9*dz**2/4/max(k_u,k_v) 
Z_range = np.arange(0,depth+dz,dz)
t_range = np.arange(0, timeTotal, dt)
Nt = int(timeTotal/dt) # total number of discrete time points
midpoint = int(Nz/2) # mid point in discrete space


## initialise 3D solution matrix
# each row&column for each different locations; each layer for each time point
global U, V
U = np.zeros((Nt+1,Nz+1, Nz+1))
V = np.zeros((Nt+1,Nz+1, Nz+1))

U[0, :, :] = np.random.rand(Nz+1,Nz+1) # initial concentration of U
V[0, :, :] = np.random.rand(Nz+1,Nz+1) # initial concentration of V


## iterate through each time point and calculate 
for i in range(1, Nt):

    #second derivative of U,V in space   
    Ut = U[i-1,0:-2,1:-1]
    Ul = U[i-1,1:-1,0:-2]
    Ub = U[i-1,2:,1:-1]
    Ur = U[i-1,1:-1,2:]
    Uc = U[i-1,1:-1,1:-1]
    depthU_2D = (Ut+Ul+Ub+Ur-4*Uc)/dz**2

    Vt = V[i-1,0:-2,1:-1]
    Vl = V[i-1,1:-1,0:-2]
    Vb = V[i-1,2:,1:-1]
    Vr = V[i-1,1:-1,2:]
    Vc = V[i-1,1:-1,1:-1]
    depthV_2D = (Vt+Vl+Vb+Vr-4*Vc)/dz**2
        
    #first derivative of U, V in time
    timeU_1D = k_u*depthU_2D + c1 - cm1*Uc + c3*Vc*Uc**2
    timeV_1D = k_v*depthV_2D + c2 - c3*Vc*Uc**2

    #calculate next time point U,V for each non-boundary point
    U[i,1:-1,1:-1] = Uc+dt*timeU_1D
    V[i,1:-1,1:-1] = Vc+dt*timeV_1D
    
    #apply boundary condition on both ends saying that no U,V going in/out
    U[i,0,:] = U[i,1,:]
    U[i,-1,:] = U[i,-2,:]
    U[i,:,0] = U[i,:,1]
    U[i,:,-1] = U[i,:,-2]

    V[i,0,:] = V[i,1,:]
    V[i,-1,:] = V[i,-2,:]
    V[i,:,0] = V[i,:,1]
    V[i,:,-1] = V[i,:,-2]

    print float(i)/float(Nt)

## display the result
zStep = int(Nz/4) # plot five locations on x-axis
xbound = np.arange(0, Nz+1, zStep)
xlabel = map(str, xbound*depth/Nz)

sampleT= int(Nt/24)
Ushow = U[0:,:,:]
#Vshow = V[0:,:,:]
Umin = Ushow.min()
Umax = Ushow.max()
#Vmin = Vshow.min()
#Vmax = Vshow.max()

fig1, ax1 = plt.subplots()
#fig2, ax2 = plt.subplots(122)

def animateU(t):
    im = plt.imshow(U[t,:,:], cmap=plt.cm.hot, vmin=Umin, vmax=Umax)
    plt.title('time = '+str(round(float(t*dt),3)))
    plt.xticks(xbound, xlabel)
    plt.yticks(xbound, xlabel) 

 
ani = animation.FuncAnimation(fig1, animateU, np.arange(1,Nt,1000),interval=1)

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

