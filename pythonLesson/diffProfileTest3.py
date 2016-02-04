import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

## This script is for testing different steady-state diffusion profiles
## we could create by having differnt sink and source efficiency

#signal saturation level
satLevel = 10

MM = 10000
## initialise parameters of heat equation
depth = 4.
timeTotal =  1000.
K = 1 # diffusion coefficient
prod = 1*MM # heat generation rate
deg = .01*MM # head absorption rate

## specify parameter for numerical solution
Nz = 200. # number of discrete step in space
dz = depth/Nz # size of discrete step in space
# array of each location in our discrete space
Z_range = np.arange(0,depth+dz,dz)
dt = 0.9*dz**2/2 # size of discrete step in time, note that dt<=dz**2/2 for stability of numerical solution
Nt = int(timeTotal/dt) # total number of discrete time points
qt1 = int(Nz/4) # location that heat gets generated
qt3 = int(3*Nz/4) # location that heat get absorbed away
qt2 = int(Nz/2) # mid point in discrete space

## initialise solution matrix
# each row for each different locations; each column for each time point
T = np.zeros((Nz+1, Nt+1)) 

print(str(dt))
## iterate through each time point and calculate 
for i in range(1, Nt):
    # second derivative of temperature in space
    depth_2D = (T[0:-2, i-1]-2*T[1:-1,i-1]+T[2:,i-1])/dz**2
    # first derivative of temperature in time
    time_1D = K*depth_2D
    # calculate next time point temperature for each non-boundary point
    T[1:-1,i] = T[1:-1,i-1]+dt*time_1D
    T[1:qt1, i] = T[1:qt1, i] + dt*prod# add heat
    T[1:Nz-1, i] = T[1:Nz-1, i] - dt*deg*T[1:Nz-1, i] # remove heat
    # apply boundary condition on both ends saying that no heat going in/out
    T[Nz,i] = T[Nz-1,i]
    T[0,i] = T[1,i]

## display the result
# we will plot heat distribution from four different time points
sampleT= int(Nt/4)
sampleTarray = sampleT*np.array([1,2,3,4])*timeTotal/Nt
sampleTarray = np.round(sampleTarray,4)
sampleTarray = map(str,sampleTarray)

#Tshow = np.flipud(T[:, sampleT:].transpose())
Tshow = np.flipud(T.transpose())


# plot location vs temperature at four different time points
plt.plot(Z_range, T[:, 1*sampleT] , color = 'r', label= 't = '+sampleTarray[0])
plt.plot(Z_range, T[:, 2*sampleT] , color = 'g', label= 't = '+sampleTarray[1])
plt.plot(Z_range, T[:, 3*sampleT] , color = 'b', label= 't = '+sampleTarray[2])
plt.plot(Z_range, T[:, 4*sampleT] , color = 'k', label= 't = '+sampleTarray[3])
plt.ylabel('temperature')
#plt.xticks([]) 
plt.legend(loc=1)
plt.xlabel('location')
plt.ylabel('temperature')
plt.axis([0, depth, -0.1, 1.1*100])

plt.show()


