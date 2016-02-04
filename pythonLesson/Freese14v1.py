import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.animation as animation

## This script simulates bacterial growth, competition, conjugation in 1D based on Freese14 Biophysical Journal paper


## initialise simulation parameter
Lsim = 1000 #1000 # number of population demes
N = 20 #100 # number of individual in each deme
Time = 5000 #100 # number of simulated generations
s = 0.01 # plasmid cost
r = 0.02 # conjugation rate
m = 0.05 # migration rate per generation.

Df = 0.5 # initial fraction of population that are donor
Rf = 0.5 # initial fraction of population that are recipient
Tf = 0 # initial fraction of population that are transconjugant

## initialise 3D matrices to store number of donor(D), recipient(R) and transconjugant(T) at different position and time points.
## ... each layer from upper to lower for donor, recipient and transconjugant
## ... each row from top to bottom for each generation time point
## ... each column for each position in 1D
DRT = np.zeros((3, Time, Lsim)) # donor

# setup numbers of D, R, T at each location for the first generation
DRT[0, 0,:] = np.round(N*Df)
DRT[1, 0,:] = np.round(N*Rf)
DRT[2, 0,:] = np.round(N*Tf)


#================ define function =============================
# this function takes a list of state transition index and returns a matrix that tells how constituent of each deme should change. 
def makeDRTchange(indList):
    # define a look up table for state transition
# ..each row from top to bottom is for donor, recipient, transconjugant change
# ..each column is for different transition paths
   stChange = np.zeros((3,7))
   stChange[:,0] = [-1,+1, 0] # r replaces d
   stChange[:,1] = [+1,-1, 0] # d replaces r
   stChange[:,2] = [ 0,+1,-1] # r replaces t
   stChange[:,3] = [ 0,-1,+1] # t replaces r
   stChange[:,4] = [+1, 0,-1] # d replaces t
   stChange[:,5] = [-1, 0,+1] # t replaces d
   stChange[:,6] = [ 0, 0, 0] # t replaces d
   
   DRTchange = np.zeros((3,Lsim))
   
   for ind in range(0,len(indList)):
     DRTchange[:,ind] = stChange[:, indList[ind]]
      
   return DRTchange

def migrate(mZone, mRate):
    # create cummulative probability table for migration
    mTable = np.array([0.5*mRate, 1-0.5*mRate])

    Nn = sum(mZone[:,0])
    #print(mZone)
    # create cummulative probability table for exchange
    eTable = np.zeros((2,3))
    eTable[0,:] = mZone[0,:]/Nn
    eTable[1,:] = mZone[0,:]/Nn+mZone[1,:]/Nn

    #get initial cell count from the center of migration 
    cellCountAll = np.copy(mZone[:,1])
    
    # go through each section of migration center. sc = 0, 1, 2 for D, R an T
    for sc in range(0,3):
        cellCount = cellCountAll[sc] # get cell number from migration center.
        
        while cellCount>0:
             destination = sum(np.random.rand() > mTable) # where the cell goes
             #type of cell from destination to return
             
             cellType = sum(np.random.rand()> eTable[:, destination])
             #print(eTable)
             #print('destination'+str(destination))
             #print('centerType'+str(sc)+'targetType'+str(cellType))
             #print(mZone)
             #print('-----------------------')
             mZone[sc,1] = mZone[sc,1]-1 # cell departs migration center
             #cell arrives at target
             mZone[sc,destination] = mZone[sc,destination]+1 
             #desnation target depart.
             mZone[cellType, destination] = mZone[cellType, destination]-1
             #destination target arrives at migration center
             mZone[cellType,1] = mZone[cellType,1] + 1
             #update eTable
             eTable[0,:] = mZone[0,:]/Nn
             eTable[1,:] = mZone[0,:]/Nn+mZone[1,:]/Nn


             #print('sc:' + str(sc) + 'cellcount:' + str(cellCount))
             cellCount = cellCount-1
    return mZone

       
#=====================================================


# iterate through each generation and calculate the number of donor, recipient, tranconjugant at each location
for i in range(1,Time):
   print('time = ' + str(i) + '/' + str(Time))
   
#  division,conj, replacement happen ~ N time in each deme in one generation
# ---------------------------------------------------------------   
   for j in range(0,N):
     # calculate fraction of donor, rec
     df = DRT[0,i-1,:]/N # donor fraction  at prev generation
     rf = DRT[1,i-1,:]/N # recipient fraction at prev  generation   
     tf = DRT[2,i-1,:]/N # recipient fraction at prev  generation

     # calculate probability of each different transition
     dtor = rf*df*(1 + 0.5*s - 0.5*r) # r replaces d
     rtod = df*rf*(1 - 0.5*s - 0.5*r) # d replaces r
     ttor = rf*tf*(1 + 0.5*s - 0.5*r) # r replaces t
     rtot = rf*df*r + rf*tf*(1 - 0.5*s + 0.5*r) # t replaces r
     ttod = tf*df # t replaces d
     dtot = df*tf # d replaces t

     # make a cummulative prop table for determining state transition 
     cumProp = np.zeros((6,Lsim))
     cumProp[0,:] = dtor
     cumProp[1,:] = dtor+rtod
     cumProp[2,:] = dtor+rtod+ttor
     cumProp[3,:] = dtor+rtod+ttor+rtot
     cumProp[4,:] = dtor+rtod+ttor+rtot+ttod
     cumProp[5,:] = dtor+rtod+ttor+rtot+ttod+dtot
   
     #generate random number array to determining state transition
     randArray = np.zeros((6,Lsim))
     randArray[:,:] = np.random.rand(Lsim)
     
     #determine state transition for each location
     ttable = randArray>cumProp
     stateIdx = ttable.sum(axis=0) # state index
   
     # DRTchanges is a 3xLsim array. Each column indicates how constituent changes in each deme. Each row indicates change in the number of d, r, t
     DRTchanges = makeDRTchange(stateIdx)
     
     # calculate and record new constituent into storage matrix
     DRT[:,i,:] = DRT[:,i-1,:]+DRTchanges



     
# ---------------------------------------------------------------
 # migration
 
   # make a temperary matrix to handle migration
   DRTiTemp =   np.copy(DRT[:,i,:])
    
   #make a random order of deme to migrate
   migrateOrder = np.array(range(0,Lsim))
   np.random.shuffle(migrateOrder)

   #iterate through the order and simulate migration
   for j in range(0,Lsim):
     # Here we get the index of migration center and indices of its left and right location. if/elif/else below are used to handle periodic boundary condition  
     migrateInd = migrateOrder[j]
     if migrateInd == 0:
        migrateZoneInd = np.array([Lsim-1, 0, 1]) #left boundary
     elif migrateInd == Lsim-1:
        migrateZoneInd = np.array([Lsim-2, Lsim-1, 0]) #right boundary
     else:
        migrateZoneInd = np.array([migrateInd-1, migrateInd, migrateInd+1])

     # migration zone is a 3x3 array. Each column indicates the number of D/R/T (top to bottom row) at the left, center and right of migration  
     migrateZone =  np.zeros((3,3))
     migrateZone[:,0] = np.copy(DRTiTemp[:, migrateZoneInd[0]])
     migrateZone[:,1] = np.copy(DRTiTemp[:, migrateZoneInd[1]])
     migrateZone[:,2] = np.copy(DRTiTemp[:, migrateZoneInd[2]])
      
     # let individual migrate locally around chosen location
     migrateZone = migrate(migrateZone, m)
     # update DRTiTemp
     DRTiTemp[:, migrateZoneInd[0]] = np.copy(migrateZone[:,0])
     DRTiTemp[:, migrateZoneInd[1]] = np.copy(migrateZone[:,1]) 
     DRTiTemp[:, migrateZoneInd[2]] = np.copy(migrateZone[:,2]) 

   #record new constituent into storage matrix
   DRT[:,i,:] = np.copy(DRTiTemp)
     

##====

#print(DRT)


DRTnew = np.zeros((Time, Lsim, 3))

DRTnew[:,:,0] = DRT[0,:,:]/N
DRTnew[:,:,2] = DRT[1,:,:]/N    
DRTnew[:,:,1] = DRT[2,:,:]/N
    
#print(DRTnew)


#plt.imshow(DRTnew)
#plt.show()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(DRTnew)

ax.set_aspect('equal')
ax.set_aspect('auto')
fig.show()
#fig.savefig('equal.png')
#ax.set_aspect('auto')
#fig.savefig('auto.png')



## display the result


              
'''
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
