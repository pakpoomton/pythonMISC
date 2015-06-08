import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

## This script shows how to solve 1D heat equation numerically.
# Here we study the progress of 1D, two species stripe pattern

#signal saturation level
satLevel = 10

## initialise parameters of heat equation
depth = 20. # in this case the width of space to generate  pattern
celSize = 1. # this is expected size of the initial band
timeTotal = 500.


D = 1 # diffusion rate of s1 and s2 signal
Ka = 0.001 # maximal switching rate of c0->c1s and c0->c2s
Km = 100 # [s2] and [s1] at half maximal swiching rate of c0->c1s and c0->c2s
A = .1 # prduction rate of s1 and s2 from from c1 and c2
Ga = .001 #maximal degradation rate s1 and s2
Gm = 1 # [s1] and [s2] at half maximal degradation rate
Kaa = .00001 # rate of conversion for c1->c1*, c2->c2*

Tau = 1.#1.99 # the delay time unit for c1s->c1 and c2s->c2 conversion


Tprior = Tau*1.1 # this is the time period before running simulation
# we wil use this time for handling delay diff eq
#**



## specify parameter for numerical solution
Nz = 400. # number of discrete step in space
dz = depth/Nz # size of discrete step in space
# array of each location in our discrete space
dt = 0.9*dz**2/2/D # size of discrete step in time, note that dt<=dz**2/4 for stability of numerical solution
Z_range = np.arange(0,depth+dz,dz)
t_range = np.arange(0, timeTotal, dt)
Nt = int(timeTotal/dt) # total number of discrete time points in simulation
NtPrior = int(Tprior/dt) # total number of discrete time prior to simulation.
cel = int(Nz*celSize/depth)



#**

## initialise solution matrix
# each row for each different locations; each column for each time point
S1 = np.zeros((Nz+1, NtPrior+Nt+1))
S2 = np.zeros((Nz+1, NtPrior+Nt+1))
C0 = np.zeros((Nz+1, NtPrior+Nt+1))
C1 = np.zeros((Nz+1, NtPrior+Nt+1))
C2 = np.zeros((Nz+1, NtPrior+Nt+1))
C1s = np.zeros((Nz+1, NtPrior+Nt+1))
C2s = np.zeros((Nz+1, NtPrior+Nt+1))

C0[cel+1:, 0:NtPrior+1] = 1 #before simulation starts, all cell at >cell is C0
C1[0:cel, :] = 1 #before simulation starts, all cell at >cell is C1



## iterate through each time point starting from the end of tprior and calculate 
for i in range(NtPrior+1,  NtPrior+Nt+1):
#for i in range(NtPrior+1,  NtPrior+5):
    
    # second derivative of each diffusible species in space
    
    depthS1_2D = (S1[0:-2, i-1]-2*S1[1:-1,i-1]+S1[2:,i-1])/dz**2
    depthS2_2D = (S2[0:-2, i-1]-2*S2[1:-1,i-1]+S2[2:,i-1])/dz**2

    S1c = S1[1:-1,i-1]
    S2c = S2[1:-1,i-1]
    # first derivative of each species concentration in time
    timeS1_1D = D*depthS1_2D
    timeS2_1D = D*depthS2_2D

    # calculate next time point cell type fraction
    C1s[:, i] =  C1s[:, i-1]+C0[:,i-1]*S2[:,i-1]*Ka/(S2[:,i-1]+Km)-Kaa*C1s[:,i-1] 
    C2s[:, i] =  C2s[:, i-1]+C0[:,i-1]*S1[:,i-1]*Ka/(S1[:,i-1]+Km)-Kaa*C2s[:,i-1] 
    C0[:,i] =   C0[:,i-1]-C0[:,i-1]*S1[:,i-1]*Ka/(S1[:,i-1]+Km) -C0[:,i-1]*S2[:,i-1]*Ka/(S2[:,i-1]+Km)
    C1[:,i] =  C1[:,i-1]+Kaa*C1s[:,i-1] 
    C2[:,i] =  C2[:,i-1]+Kaa*C2s[:,i-1]

    #print [C2s[30, i], C2[30, i]]
    #print C0[30,i-1]
    #print C0[30,i-1]*S2[30,i-1]*Ka/(S2[30,i-1]+Km)
    #print S1[30,i-1]*Ka/(S1[30,i-1]+Km)
    #print S1[30,i-1]
    
    # calculate next time point concentrations for each non-boundary point
    S1[1:-1,i] = S1c + dt*timeS1_1D 
    S2[1:-1,i] = S2c + dt*timeS2_1D

    S1[1:-1,i] = S1[1:-1,i] + A*C1[1:-1,i] - Ga*S1[1:-1,i]
    S2[1:-1,i] = S2[1:-1,i] + A*C2[1:-1,i] - Ga*S2[1:-1,i]
    

    # apply boundary condition on both ends saying that nothing going in/out
    S1[Nz,i] = S1[Nz-1,i]
    S2[Nz,i] = S2[Nz-1,i]
    S1[0,i] = S1[1,i]
    S2[0,i] = S2[1,i]

## display the result
# we will plot heat distribution from four different time points
sampleT= int(Nt/4)
sampleTs = sampleT*np.array([0.5, 2, 3 ,4])
sampleTarray = sampleTs*timeTotal/Nt
sampleTarray = np.round(sampleTarray,2)
strTarray = map(str,sampleTarray)

S1show = np.flipud(S1[:, NtPrior:].transpose())


# plot location vs S1 at four different time points
plt.subplot(411)
plt.plot(C0[:, -1],color='k')
plt.plot(C1[:, -1],color='r')
plt.plot(C2[:, -1],color='b')
plt.plot(C1s[:, -1],color='r', linestyle='--')
plt.plot(C2s[:, -1],color='b', linestyle='--')

plt.ylabel('fraction')
plt.axis([0, Nz+1, -0.1, 1.1])
plt.xticks([]) 


# plot location vs S1 at four different time points
plt.subplot(412)
plt.plot(Z_range, S1[:, NtPrior+1*sampleT],color='r',label= 't = '+strTarray[0])
plt.plot(Z_range, S1[:, NtPrior+2*sampleT],color='g',label= 't = '+strTarray[1])
plt.plot(Z_range, S1[:,NtPrior+3*sampleT],color='b',label= 't = '+strTarray[2])
plt.plot(Z_range, S1[:,NtPrior+4*sampleT],color= 'k',label= 't = '+strTarray[3])
plt.ylabel('concentration')
plt.xticks([]) 
plt.legend(loc=1)

# plot location vs S1 at four different time points
plt.subplot(413)
plt.plot(Z_range, S1[:, NtPrior+1*sampleT],color='r')
plt.plot(Z_range, S1[:, NtPrior+2*sampleT],color='g')
plt.plot(Z_range, S1[:,NtPrior+3*sampleT],color='b')
plt.plot(Z_range, S1[:,NtPrior+4*sampleT],color= 'k')
plt.ylabel('concentration')
plt.xticks([]) 
plt.axis([0, depth, 0, satLevel])


# make heatmap of S1 in location and time
plt.subplot(414)
plt.imshow(S1show, cmap=plt.cm.hot, interpolation='nearest', origin='lower', vmin=0, vmax=satLevel)
plt.axis('auto')
plt.axis([0, Nz+1, 0, Nt+1])
zStep = int(Nz/5) # plot five locations on x-axis
xbound = np.arange(0, Nz+1, zStep)
xlabel = map(str, xbound*depth/Nz)
ybound = [0, Nt]
timeStartDisplay =  round(timeTotal*sampleT/Nt, 2)
timeStartDisplay =  0
ylabel = map(str, [timeTotal, timeStartDisplay])
plt.xticks(xbound, xlabel) 
plt.yticks(ybound, ylabel)
plt.xlabel('location')
plt.ylabel('time')
# setup color bar
S1min = S1show.min()
S1max = S1show.max()
S1step = (S1max-S1min)/4
colorRange = np.arange(S1min, S1max+S1step, S1step)
colorRange = np.round(colorRange, 2)
plt.colorbar(orientation='horizontal', ticks = colorRange)



plt.show()

#print C1s[:, -10]
#print C2s[:, -10]

#print C1[:, -10]
#print C2[:, -10]
