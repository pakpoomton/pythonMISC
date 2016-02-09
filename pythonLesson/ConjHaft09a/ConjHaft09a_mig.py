# plasmid conjugation  modeling based on Haft et al 2009 "Competition favour redced cost of plasmids to host bacteria
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import time

startTime = time.time() # time when simulation starts

# parameter set based on Table 1 of Haft et al 2009
N_0 = 10**6.5 # initial density of plasmid-free cells (CFU/ml)
P1a_0 = 0 # initial density of cells w/ transitory-derepressed Pfin+ (CFU/ml)
P1b_0 = 10**4.5 # initial density of cells w/ Pfin+ (CFU/ml)
P2_0 = 10**4.5 # initial density of cells w/ Pfin- (CFU/ml)
C_0 = 200 # initial conc of resource, C (ug/ml)
rN = 1.459 # growth rate of N (/h)
r1a = 1.230 # growth rate of P1a (/h)
r1b = 1.405 #  growth rate of P1b (/h)
r2 =  1.230 # growth rate of P2 (/h)
y1a = 3.8e-9 # conjugation rate of P1a donors (ml/cell/h)
y1b = 4.4e-12 # conjugation rate of P1b donors (ml/cell/h)
y2 = 3.8e-9 # conjugation rate of P2 donors (ml/cell/h)
s1 = 1e-4 # segregation rate for Pfin+
s2 = 1e-4 # segregation rate for Pfin-
f1 = 0.1 # fin repression rate for Pfin+
Km = 0.2 # Monod constant (ug/ml)
Y = 8e-8 # yield coefficient (ug/CFU)

N_mig = 10**6.4 # plasmid free cell migration each dilution cycle

Cyc = 100 # number of growth cycles (i.e. 24 hr dilution cycle)
# ** try 500 cycles to see oscillation!
dFold = 1000. # dilution fold for each growth cycle

# setup an initial condition vector
initState = [N_0, P1a_0, P1b_0, P2_0, C_0]

# setup time vector (h) for each growth cycle
t = np.linspace(0, 24., 1000) # time grid

solnALL = np.copy([initState]) # array to store all time course solution
solnEndDay = np.copy([initState]) # array to store end of day condition
tALL = np.linspace(t[0], t[-1]*Cyc, len(t)*Cyc) # all time array

# define growth equation (growth depend on available nutrient)
def grow(nutrient):
        return nutrient/(Km+nutrient)

# define conj rate scaling equation
def conj(nutrient):
        if (nutrient > 0.1*Km):
                return 1
        else :
                return 0.001

# solve the system dy/dt = f(y, t)
def dState(state,time):
        N = state[0]
        P1a = state[1]
        P1b = state[2]
        P2 = state[3]
        C = state[4]

        # ODEs from Haft et al 2009
        dState0 = (rN*N + r1a*s1*P1a + r1b*s1*P1b + r2*s2*P2)*grow(C) - \
                  (y1a*P1a + y1b*P1b + y2*P2)*N*conj(C)
        dState1 = (1-s1)*(r1a*P1a)*grow(C) + \
                  (y1a*P1a*N)*conj(C) + (y1b*P1b*N)*conj(C) - (f1*P1a)
        dState2 = (1-s1)*(r1b*P1b)*grow(C) + (f1*P1a)
        dState3 = (1-s2)*(r2*P2)*grow(C) + (y2*P2*N)*conj(C)
        dState4 = -Y*grow(C)*(rN*N + r1a*P1a + r1b*P1b + r2*P2)
        
        return [dState0, dState1, dState2, dState3, dState4]


startState = np.copy(initState)

for k in range(0, Cyc):
        soln = odeint(dState, startState, t) # solve ODEs for this cycle
        solnALL = np.concatenate((solnALL, soln), axis=0) # store result
        
        # get ODEs result at the end of the day and store 
        solnEnd = np.array([soln[-1,:]]) 
        solnEndDay = np.concatenate((solnEndDay, solnEnd), axis=0)

        # no. of each cell type and resource level at the end of the cycle 
        N_final =  soln[-1,0]
        P1a_final = soln[-1,1]
        P1b_final = soln[-1,2]
        P2_final = soln[-1,3]
        C_final = soln[-1,4]
        
        # dilute the sample and regrow in fresh media at the end of the cycle
        startState = [N_final/dFold + N_mig*(dFold-1)/dFold, \
                      P1a_final/dFold, P1b_final/dFold, \
                      P2_final/dFold, C_final/dFold + \
                      C_0*(dFold-1)/dFold]
          
########


N_array = solnALL[1:,0]
P1_array = solnALL[1:,1] + solnALL[1:,2]
P2_array = solnALL[1:,3]
C_array = solnALL[1:,4]/Y # scale food to unit ~ cell

Nd_array = solnEndDay[:,0]
P1d_array = solnEndDay[:,1] + solnEndDay[:,2]
P2d_array = solnEndDay[:,3]
Cd_array = solnEndDay[:,4]/Y # scale food to unit ~ cell

N_array = np.log10(N_array)
P1_array = np.log10(P1_array)
P2_array = np.log10(P2_array)
C_array = np.log10(C_array)


Nd_array = np.log10(Nd_array)
P1d_array = np.log10(P1d_array)
P2d_array = np.log10(P2d_array)
Cd_array = np.log10(Cd_array)

########
dayArray = np.array(range(0,Cyc+1)) # end of cycle array index (for plotting)

print time.time()-startTime, "seconds wall time"

# plot results
plt.subplot(211)
plt.plot(tALL, N_array, label='Free')
plt.plot(tALL, P1_array, label='fin+')
plt.plot(tALL, P2_array, label='fin-')
#plt.plot(tALL, C_array, label='food')
plt.xlabel('time(hr)')
plt.ylabel('log10(CFU/mL)')
plt.title('Plasmid Infection Dynamics')
plt.axis([0, 2500, 0, 10])
plt.legend(loc=4)


plt.subplot(212)
plt.plot(dayArray, Nd_array, "D-", label='Free')
plt.plot(dayArray, P1d_array, "s-", label='fin+')
plt.plot(dayArray, P2d_array, "^-", label='fin-')
#plt.plot(dayArray, Cd_array, label='food')
plt.xlabel('time (day)')
plt.ylabel('log10(CFU/mL)')
plt.axis([0, 100, 0, 10])
plt.legend(loc=4)

plt.show()


