# plasmid conjugation  modeling
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
#plt.ion() # this allow interactive mode so we can input through command line
          # and see changes in the figure right away without rerunning it.

#define parameters for conjugation and growth
e = 1 # resource required per uninfected cell
e_plus = 1 # resource required per infected cell
psi = 1
psi_plus = 1
gamma = .01 # conj rate
tau = 0 # segregation loss rate

gSc = 1 # fitness scaling for uninfected cells
gSc_plus = .7 # fitness scaling for infected cells

P = 1 # max growth rate
Q = 1 # resource level at half growth rate

Cyc = 100 # number of growth cycles
dFold = 100. # dilution fold for each growth cycle

# initial conditions
r0 = 100        # initial nutrient level
n_plus0 = 0.01    # initial cell bearing conjugation plasmid
n0 = 1         # initial cell without conjugation plasmid
y0 = [r0, n_plus0, n0] # initial condition vector
t = np.linspace(0, 5., 1000) # time grid
solnALL = np.copy([y0]) # array to store all time course solution
solnEndDay = np.copy([y0]) # array to store condition at the end of the day
tALL = np.linspace(t[0], t[-1]*Cyc, len(t)*Cyc) 

# define growth equation
def g(nutrient):
        return P*nutrient/(Q+nutrient)

# solve the system dy/dt = f(y, t)
def f(y,t):
        r = y[0]
        n_plus = y[1]
        n = y[2]
        # the model equation adapted from Stward &  Levin 1977
        f0 = -e_plus*g(r)*n_plus - e*g(r)*n
        f1 = n_plus*g(r)*gSc_plus + gamma*n_plus*n - tau*n_plus
        f2 = n*g(r)*gSc - gamma*n_plus*n + tau*n_plus
        return [f0, f1, f2]


yStart = np.copy(y0)

for k in range(0, Cyc):
  # solve the DEs
  soln = odeint(f, yStart, t)
  solnALL = np.concatenate((solnALL, soln), axis=0)
  solnEnd = np.array([soln[-1,:]])
  solnEndDay = np.concatenate((solnEndDay, solnEnd), axis=0)
  r_final = soln[-1,0]
  n_plus_final = soln[-1,1]
  n_final = soln[-1,2]
  yStart = [r_final/dFold + r0*(dFold-1)/dFold, n_plus_final/dFold, n_final/dFold]
  
r_array = solnALL[1:, 0]
n_plus_array = solnALL[1:, 1]
n_array = solnALL[1:, 2]

r_EndDay = solnEndDay[:, 0]
n_plus_EndDay = solnEndDay[:, 1]
n_EndDay = solnEndDay[:, 2]


#print(r_array.shape)
#print(tALL.shape)

dayArray = range(0,Cyc+1)

# plot results
plt.subplot(211)
plt.plot(tALL, r_array, label='Nutrient')
plt.plot(tALL, n_plus_array, label='Infected')
plt.plot(tALL, n_array, label='Free')
plt.xlabel('time')
plt.ylabel('Population')
plt.title('Plasmid Infection Dynamics')
plt.legend(loc=1)

plt.subplot(212)
plt.plot(dayArray, r_EndDay, label='Nutrient')
plt.plot(dayArray, n_plus_EndDay, label='Infected')
plt.plot(dayArray, n_EndDay, label='Free')
plt.xlabel('time (cycle)')
plt.ylabel('Population')
plt.legend(loc=1)

plt.show()


