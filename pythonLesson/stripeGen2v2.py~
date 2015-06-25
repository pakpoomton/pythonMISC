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
timeTotal = 10.


D = 1 # diffusion rate of s1 and s2 signal
Ka = 0.1 # maximal switching rate of c0->c1s and c0->c2s
Km = 10 # [s2] and [s1] at half maximal swiching rate of c0->c1s and c0->c2s
A = .1 # prduction rate of s1 and s2 from from c1 and c2
Ga = .001 #maximal degradation rate s1 and s2
Gm = 1 # [s1] and [s2] at half maximal degradation rate

Tau = 1# the delay time unit for c1s->c1 and c2s->c2 conversion






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
C1show = np.flipud(C1.transpose())
C2show = np.flipud(C2.transpose())

C1sSum = np.sum(C1s,axis=1)
C2sSum = np.sum(C2s,axis=1)

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

# plot location vs S1 at four different time points
plt.subplot(411)
#plt.imshow(C1show, cmap=plt.cm.gray, vmin=0, vmax=1)
plt.imshow(C1C2)
plt.axis('auto')
plt.axis([0, Nz+1, 0, Nt+1])
plt.ylabel('time')

'''
plt.plot(C0[:, -1],color='k')
plt.plot(C1[:, -1],color='r')
plt.plot(C2[:, -1],color='b')
plt.plot(C1sSum,color='r', linestyle='--')
plt.plot(C2sSum,color='b', linestyle='--')

plt.ylabel('fraction')
plt.axis([0, Nz+1, -0.1, 1.1])
plt.xticks([]) 
'''

# plot location vs S1 at four different time points
plt.subplot(412)
plt.plot(Z_range, S1[:, 1*sampleT],color='r',label= 't = '+strTarray[0])
plt.plot(Z_range, S1[:, 2*sampleT],color='g',label= 't = '+strTarray[1])
plt.plot(Z_range, S1[:, 3*sampleT],color='b',label= 't = '+strTarray[2])
plt.plot(Z_range, S1[:, 4*sampleT],color= 'k',label= 't = '+strTarray[3])
plt.ylabel('concentration')
plt.xticks([]) 
plt.legend(loc=1)

# plot location vs S1 at four different time points
plt.subplot(413)
plt.plot(Z_range, S1[:, 1*sampleT],color='r')
plt.plot(Z_range, S1[:, 2*sampleT],color='g')
plt.plot(Z_range, S1[:, 3*sampleT],color='b')
plt.plot(Z_range, S1[:, 4*sampleT],color= 'k')
plt.ylabel('concentration')
plt.xticks([]) 
plt.axis([0, depth, 0, satLevel])


# make heatmap of S1 in location and time
plt.subplot(414)
plt.imshow(S1show, cmap=plt.cm.hot, interpolation='nearest', origin='lower', vmin=0, vmax=satLevel)
plt.axis('auto')
plt.axis([0, Nz+1, 0, Nt+1])

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

'''
print C1sSum
print C2sSum

print len(C1sSum)
print len(C2sSum)
print len(C0[:,-1])
print C0[:, -1] + C1[:, -1] + C2[:, -1] + C1sSum + C2sSum
'''
#print np.sum(C1s,axis=1)
#print np.sum(C2s,axis=1)
#print len(np.sum(C1s,axis=1))
#print Nz

#print C0[50,Nt]
#print C1[50,Nt]
#print C2[50,Nt]
#print sum(C1s[50,:])
#print sum(C2s[50,:])

#print C0[:, 49]

#print C1s[:, -10]
#print C2s[:, -10]

#print C1[:, -10]
#print C2[:, -10]
