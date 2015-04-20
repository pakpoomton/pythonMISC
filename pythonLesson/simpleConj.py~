# plasmid conjugation  modeling
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
plt.ion() # this allow interactive mode so we can input through command line
          # and see changes in the figure right away without rerunning it.

#define parameters for conjugation and growth
e = 1
e_plus = 1
psi = 1
psi_plus = 1
gamma = .1
tau = 0

P = 1
Q = 1


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
        f1 = n_plus*g(r) + gamma*n_plus*n - tau*n_plus
        f2 = n*g(r) - gamma*n_plus*n + tau*n_plus
        return [f0, f1, f2]

# initial conditions
r0 = 100        # initial nutrient level
n_plus0 = 0.01    # initial cell bearing conjugation plasmid
n0 = 1         # initial cell without conjugation plasmid
y0 = [r0, n_plus0, n0] # initial condition vector
t = np.linspace(0, 5., 1000) # time grid


# solve the DEs
soln = odeint(f, y0, t)
r_array = soln[:, 0]
n_plus_array = soln[:, 1]
n_array = soln[:, 2]

# plot results
plt.figure()
plt.plot(t, r_array, label='Nutrient')
plt.plot(t, n_plus_array, label='Infected')
plt.plot(t, n_array, label='Free')
plt.xlabel('time')
plt.ylabel('Population')
plt.title('Plasmid Infection Dynamics')
plt.legend(loc=0)


# Two subplots, the axes array is 1-d
'''
f, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(t, r_array)
axarr[0].set_title('Plasmid Infection Dynamics')
#axarr[0].ylabel('Nutrient')
axarr[1].scatter(t, n_plus_array, label='Infected')
axarr[1].scatter(t, n_array, label='Free')
#axarr[1].ylabel('Population')
'''
