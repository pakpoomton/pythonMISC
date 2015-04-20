import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt



## set-up simulation
# define simulation time
start = 0 # starting time point
end = 2   # final time point
numsteps = 1000 # time step
time = np.linspace(start,end,numsteps)

# set-up initial condition (initial population)
y0=np.array([10])


# define exponential growth function 
r  = 2 # intrinsic growth rate
K = 100 # carrying capacity

def deriv(y,t):
    yprime = np.array([r*y[0]])
    return yprime

# define logistic exponential growth function 
def deriv2(y,t):
    yprime = np.array([r*y[0]*(1-y[0]/K)])  
    return yprime

# solve ODEs
y=odeint(deriv, y0, time) # solve ODE for exponential growth
y2=odeint(deriv2, y0, time) # solve ODE for logistic growth

# generate plots
plt.plot(time,y[:], color = 'r', label='exponent')
plt.plot(time,y2[:], color = 'b', label = "logistic")
plt.title('Population Growth')
plt.legend(loc = 2)
plt.xlabel('Time')
plt.ylabel('Population')
plt.show()
