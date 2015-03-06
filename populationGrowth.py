import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

r  = 2 # intrinsic growth rate
K = 100 # carrying capacity

# define exponential growth function 
def deriv(y,t):
    yprime = np.array([intrinsicGrowthRate*y[0]])
    return yprime

# define logistic exponential growth function 
def deriv2(y,t):
    yprime = np.array([intrinsicGrowthRate*y[0]*(1-y[0]/K)])
    #yprime = np.array([intrinsicGrowthRate*y[0]*(1)])   
    return yprime


# define simulation time
start = 0
end = 1
numsteps = 1000
time = np.linspace(start,end,numsteps)
y0=np.array([10])

y=integrate.odeint(deriv, y0, time)
line_up = plt.plot(time,y[:], color = 'r', label='exponent')

y2=integrate.odeint(deriv2, y0, time)
line_down = plt.plot(time,y2[:], color = 'b', label = "logistic")
plt.title('Population Growth')
plt.legend(loc = 2)
plt.xlabel('Time')
plt.ylabel('Population')
plt.show()
