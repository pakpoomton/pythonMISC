import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import time

start = time.time()
## This script shows how to solve 1D heat equation numerically.
# Here we study the progress of 1D, two species stripe pattern



## initialise parameters of heat equation
depth = 30. # in this case the width of space to generate  pattern
celSize = 1. # this is expected size of the initial band
timeTotal = 15.


D = 10 # diffusion rate of s1 and s2 signal
Ka = 1 # maximal switching rate of c0->c1s and c0->c2s
Km = 0.1 # [s2] and [s1] at half maximal swiching rate of c0->c1s and c0->c2s
A = .1 # production rate of s1 and s2 from from c1 and c2
Ga = .1 #maximal degradation rate s1 and s2
Gm = 1 # [s1] and [s2] at half maximal degradation rate

Tau = 1# the delay time unit for c1s->c1 and c2s->c2 conversion

#signal saturation level for making heatmap
satLevel = Km*2


## specify parameter for numerical solution
Nz = 400. # number of discrete step in space
dz = depth/Nz # size of discrete step in space
# array of each location in our discrete space
dt = 0.9*dz**2/2/D # size of discrete step in time, note that dt<=dz**2/4 for stability of numerical solution
Z_range = np.arange(0,depth+dz,dz)
t_range = np.arange(0, timeTotal, dt)
Nt = int(timeTotal/dt) # total number of discrete time points in simulation
NTau = int(Tau/dt) # total number of discrete time for production delay.
cel = int(Nz*celSize/depth)



## initialise solution matrix
# each row for each different locations; each column for each time point
S1 = np.zeros((Nz+1, Nt+1))
S2 = np.zeros((Nz+1, Nt+1))
C0 = np.zeros((Nz+1, Nt+1))
C1 = np.zeros((Nz+1, Nt+1))
C2 = np.zeros((Nz+1, Nt+1))
C1s = np.zeros((Nz+1, NTau+1))
C2s = np.zeros((Nz+1, NTau+1))

C0[cel+1:, 0] = 1 #before simulation starts, all cell at >cell is C0
C1[0:cel+1, :] = 1 #before simulation starts, all cell at >cell is C1



## iterate through each time point starting from the end of tprior and calculate 
for i in range(1,  Nt+1):
#for i in range(1,  50):
    
    # second derivative of each diffusible species in space
    
    depthS1_2D = (S1[0:-2, i-1]-2*S1[1:-1,i-1]+S1[2:,i-1])/dz**2
    depthS2_2D = (S2[0:-2, i-1]-2*S2[1:-1,i-1]+S2[2:,i-1])/dz**2

    S1c = S1[1:-1,i-1]
    S2c = S2[1:-1,i-1]
    # first derivative of each species concentration in time
    timeS1_1D = D*depthS1_2D
    timeS2_1D = D*depthS2_2D

    
    # calculate next time point cell type fraction
    outC1s = np.copy(C1s[:, NTau])
    oldC1s = np.copy(C1s[:, 0:NTau])
    inC1s = C0[:,i-1]*S2[:,i-1]*Ka/(S2[:,i-1]+Km)
    C1s[:,0] = np.copy(inC1s)
    C1s[:,1:NTau+1] = np.copy(oldC1s)

    

    outC2s = np.copy(C2s[:, NTau])
    #print C2s[50,NTau]
    oldC2s = np.copy(C2s[:, 0:NTau])
    #print C2s[50, 0:NTau]
    inC2s = C0[:,i-1]*S1[:,i-1]*Ka/(S1[:,i-1]+Km)
    #print inC2s[50]
    C2s[:,0] = np.copy(inC2s)
    C2s[:,1:NTau+1] = np.copy(oldC2s)
    #print C2s[50,:]

    #print "---"
    
    C0[:,i] =   C0[:,i-1]-inC1s-inC2s
    C1[:,i] =  C1[:,i-1]+outC1s 
    C2[:,i] =  C2[:,i-1]+outC2s

    '''
    ind = int(4./dz)
    print C0[ind,i]
    print C1[ind,i]
    print C2[ind,i]
    print sum(C1s[ind,:])
    print sum(C2s[ind,:])
    print "---"
    '''
    
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

S1show = np.flipud(S1.transpose())
S2show = np.flipud(S2.transpose())
C1show = np.flipud(C1.transpose())
C2show = np.flipud(C2.transpose())

C1sSum = np.around(np.sum(C1s,axis=1), decimals = 3)
C2sSum = np.around(np.sum(C2s,axis=1), decimals = 3)

#make image from C1 and C2 level
pp = C1show.shape
C1C2 = np.zeros((pp[0], pp[1], 3))
C1C2[:, :, 0] = np.copy(C1show)
C1C2[:, :, 2] = np.copy(C2show)

zStep = int(Nz/5) # plot five locations on x-axis
xbound = np.arange(0, Nz+1, zStep)
xlabel = map(str, xbound*depth/Nz)
ybound = [0, Nt]
timeStartDisplay =  round(timeTotal*sampleT/Nt, 2)
timeStartDisplay =  0
ylabel = map(str, [timeTotal, timeStartDisplay])


C0final =  np.around(C0[:, -1], decimals = 3)
C1final = np.around(C1[:, -1], decimals = 3)
C2final =  np.around(C2[:, -1], decimals = 3)

# Generate plots
plt.subplot(411) #fraction of different cell types at the final time point
plt.plot(C0final,color='k')
plt.plot(C1final,color='r')
plt.plot(C2final,color='b')
plt.plot(C1sSum,color='r', linestyle='--')
plt.plot(C2sSum,color='b', linestyle='--')
plt.ylabel('fraction')
plt.axis([0, Nz+1, -0.1, 1.1])
plt.xticks([])



plt.subplot(412) #[S1] w.r.t time and location

plt.imshow(S1show, cmap=plt.cm.gray, interpolation='nearest', origin='lower', vmin=0, vmax=satLevel)
plt.axis('auto')
plt.axis([0, Nz+1, 0, Nt+1])
plt.xticks([])
plt.yticks(ybound, ylabel)
plt.ylabel('time')

plt.subplot(413) #[S2] w.r.t time and location

plt.imshow(S2show, cmap=plt.cm.gray, interpolation='nearest', origin='lower', vmin=0, vmax=satLevel)
plt.axis('auto')
plt.axis([0, Nz+1, 0, Nt+1])
plt.xticks([])
plt.yticks(ybound, ylabel)
plt.ylabel('time')

# make heatmap of C1 and C2 w.r.t location and time
plt.subplot(414)

plt.imshow(C1C2, vmin =0, vmax =1)
plt.axis('auto')
plt.axis([0, Nz+1, 0, Nt+1])
plt.xticks(xbound, xlabel) 
plt.yticks(ybound, ylabel)
plt.xlabel('location')
plt.ylabel('time')

end = time.time()


print('runtime = ' + str(end-start))
plt.show()


