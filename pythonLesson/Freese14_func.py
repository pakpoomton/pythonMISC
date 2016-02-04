import numpy as np
import scipy as sp

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

    N = sum(mZone[0,:])
    # create cummulative probability table for exchange
    eTable = np.zeros((2,3))
    eTable[0,:] = mZone[0,:]/N
    eTable[1,:] = mZone[0,:]+mZone[1,:]/N

    #get initial cell count from the center of migration 
    cellCountAll = np.copy(mZone[:,1])
    
    # go through each section of migration center. sc = 0, 1, 2 for D, R an T
    for sc in range(0,3):
        cellCount = cellCountAll[sc] # get cell number from migration center.
        
        while cellCount>0:
             destination = sum(np.random.rand() > mTable) # where the cell goes
             #type of cell from destination to return
             cellType = sum(np.random.rand()> eTable[:, destination])  
             mZone[sc,1] = mZone[sc,1]-1 # cell departs migration center
             #cell arrives at target
             mZone[sc,destination] = mZone[sc,destination]+1 
             #desnation target depart.
             mZone[cellType, destination] = mZone[cellType, destination]-1
             #destination target arrives at migration center
             mZone[cellType,1] = mZone[cellType,1] + 1
             #update eTable
             eTable[0,:] = mZone[0,:]/N
             eTable[1,:] = mZone[0,:]+mZone[1,:]/N

             print('sc:' + str(sc) + 'cellcount:' + str(cellCount))
             cellCount = cellCount-1
    return mZone
       
