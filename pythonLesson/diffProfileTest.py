import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

## This script is for testing different steady-state diffusion profiles
## we could create by having differnt sink and source efficiency

#signal saturation level
satLevel = 10

## initialise parameters of heat equation
depth = 4.
timeTotal = 1.
K = 1 # diffusion coefficient
prod = .1 # heat generation rate
deg = 0.001#.01 # head absorption rate

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


## iterate through each time point and calculate 
for i in range(1, Nt):
    # second derivative of temperature in space
    depth_2D = (T[0:-2, i-1]-2*T[1:-1,i-1]+T[2:,i-1])/dz**2
    # first derivative of temperature in time
    time_1D = K*depth_2D
    # calculate next time point temperature for each non-boundary point
    T[1:-1,i] = T[1:-1,i-1]+dt*time_1D
    T[1:qt1, i] = T[1:qt1, i] + prod# add heat
    T[1:Nz-1, i] = T[1:Nz-1, i] - deg*T[1:Nz-1, i] # remove heat
    # apply boundary condition on both ends saying that no heat going in/out
    T[Nz,i] = T[Nz-1,i]
    T[0,i] = T[1,i]

## display the result
# we will plot heat distribution from four different time points
sampleT= int(Nt/4)
sampleTarray = sampleT*np.array([1,2,3,4])*timeTotal/Nt
sampleTarray = np.round(sampleTarray,2)
sampleTarray = map(str,sampleTarray)

#Tshow = np.flipud(T[:, sampleT:].transpose())
Tshow = np.flipud(T.transpose())


# plot location vs temperature at four different time points
plt.subplot(311)
plt.plot(Z_range, T[:, 1*sampleT] , color = 'r', label= 't = '+sampleTarray[0])
plt.plot(Z_range, T[:, 2*sampleT] , color = 'g', label= 't = '+sampleTarray[1])
plt.plot(Z_range, T[:, 3*sampleT] , color = 'b', label= 't = '+sampleTarray[2])
plt.plot(Z_range, T[:, 4*sampleT] , color = 'k', label= 't = '+sampleTarray[3])
plt.ylabel('temperature')
plt.xticks([]) 
plt.legend(loc=1)


# plot location vs temperature at four different time points
plt.subplot(312)
plt.plot(Z_range, T[:, 1*sampleT] , color = 'r')
plt.plot(Z_range, T[:, 2*sampleT] , color = 'g')
plt.plot(Z_range, T[:, 3*sampleT] , color = 'b')
plt.plot(Z_range, T[:, 4*sampleT] , color = 'k')
plt.ylabel('temperature')
plt.xticks([]) 
plt.axis([0, depth, 0, satLevel])



# make heatmap of temperature in location and time
plt.subplot(313)
plt.imshow(Tshow, cmap=plt.cm.hot, interpolation='nearest', origin='lower', vmin = 0, vmax = satLevel)
#plt.imshow(Tshow, cmap=plt.cm.hot, interpolation='nearest', origin='lower')

plt.axis('auto')
plt.axis([0, Nz+1, 0, Nt+1])
zStep = int(Nz/4) # plot five locations on x-axis
xbound = np.arange(0, Nz+1, zStep)
xlabel = map(str, xbound*depth/Nz)
ybound = [0, Nt-sampleT]
tSD =  round(timeTotal*sampleT/Nt, 2)
Tdisplay = [0, tSD, 2*tSD, 3*tSD, 4*tSD]
Ttickloc = [4*sampleT, 3*sampleT, 2*sampleT, 1*sampleT, 0]
TdisplayTxt = map(str, Tdisplay)
plt.xticks(xbound, xlabel) 
plt.yticks(Ttickloc, TdisplayTxt)
plt.xlabel('location')
plt.ylabel('time')
# setup color bar
Tmin = Tshow.min()
Tmax = Tshow.max()
Tstep = (Tmax-Tmin)/4
colorRange = np.arange(Tmin, Tmax+Tstep, Tstep)
colorRange = np.round(colorRange, 2)
plt.colorbar(orientation='horizontal', ticks = colorRange)

plt.show()


