# plasmid conjugation  modeling
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
#plt.ion() # this allow interactive mode so we can input through command line
          # and see changes in the figure right away without rerunning it.

#define parameters for conjugation and growth
e = 1 # nutrient take up rate for uninfected cell
e_plus = 1 # nutrient take up rate for infected cell 
psi = 1 # replication rate for uninfected cell
psi_plus = 1 #replication rate for infected cell
#gamma = .1 #conjugation rate
tau = 0 # plasmid loss rate

P = 1 # maximal replication rate
Q = 1 # dissociation constant for untrient take up for cell growth


AHL_prod = 1 # signaling molecule production rate
gamma_max = 0.1 # maximal conjugation rate (when AHL level is saturated)
AHL_diss = 5 # dissociation eq constant for AHL to conj machinery 
AHL_n = 100 # AHL cooperativity

AHL_degMax = 0 # maximal degradation rate of AHL by AiiA
AHL_degBind = 1 # binding dissociation const for AHL to deg machinery

# define growth equation
def g(nutrient):
        return P*nutrient/(Q+nutrient)

# define conjugation rate equation
def gamma(signal):
        return gamma_max*AHL_diss**AHL_n/(AHL_diss**AHL_n + signal**AHL_n)

# define degradation process of AHL
def deg(signal):
        return signal*AHL_degMax/(signal + AHL_degBind)

# solve the system dy/dt = f(y, t)
def f(y,t):
        r = y[0]
        n_plus = y[1]
        n = y[2]
        ahl = y[3]
        # the model equation adapted from Stward &  Levin 1977
        f0 = -e_plus*g(r)*n_plus - e*g(r)*n
        f1 = n_plus*g(r) + gamma(ahl)*n_plus*n - tau*n_plus
        f2 = n*g(r) - gamma(ahl)*n_plus*n + tau*n_plus
        f3 = n_plus*AHL_prod - n*deg(ahl)
        return [f0, f1, f2, f3]

# initial conditions
r0 = 100        # initial nutrient level
n_plus0 = 0.01    # initial cell bearing conjugation plasmid
n0 = 1         # initial cell without conjugation plasmid
ahl0 = 0

#r0 = 0        # initial nutrient level
#n_plus0 = 40    # initial cell bearing conjugation plasmid
#n0 = 20         # initial cell without conjugation plasmid
#ahl0 = 0



y0 = [r0, n_plus0, n0, ahl0] # initial condition vector
t = np.linspace(0, 10., 1000) # time grid


# solve the DEs
soln = odeint(f, y0, t)
r_array = soln[:, 0]
n_plus_array = soln[:, 1]
n_array = soln[:, 2]
AHL_array = soln[:, 3]

# plot results
'''
plt.figure()
plt.plot(t, r_array, label='Nutrient')
plt.plot(t, n_plus_array, label='Infected')
plt.plot(t, n_array, label='Free')
plt.xlabel('time')
plt.ylabel('Population')
plt.title('Plasmid Infection Dynamics')
plt.legend(loc=0)
'''


# Two subplots, the axes array is 1-d
f, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(t, AHL_array, color = 'k', marker = '.')
axarr[0].set_title('Plasmid Infection Dynamics')
axarr[0].set_ylabel('AHL')
#axarr[0].axis([min(t), max(t), -0.1, 5*AHL_diss])
axarr[1].scatter(t, r_array, label='Nutrient', color = 'b', marker = '.')
axarr[1].scatter(t, n_plus_array, label='Infected', color = 'g', marker = '.')
axarr[1].scatter(t, n_array, label='Free', color = 'r', marker = '.')
axarr[1].axis([min(t), max(t), -1, r0+n_plus0+n0])
plt.legend(loc=2)
axarr[1].set_ylabel('Population')
axarr[1].set_xlabel('Time')
