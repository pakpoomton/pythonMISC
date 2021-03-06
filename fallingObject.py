from numpy import *
import matplotlib.pyplot as plt

def fallingObject(initialVelocity):
    t_max = 10 #total simulation time
    t_step = 0.1 #discrete step size in simulation
    g = -9.81 #gravity
    n_steps = int(t_max/t_step)
    t = zeros(n_steps)
    v = zeros(n_steps)

    #set initial condition
    t[0] = 0 # initial time point
    v[0] = initialVelocity; # initial velocity

    # solving for velocity on each time point using Euler method
    for i in range(1, n_steps):
        t[i] = t[i-1] + t_step
        v[i] = v[i-1] + t_step*g

    print("time = " + str(t))
    print("velocity = " + str(v))

    plt.plot(t, v, linewidth = 2)
    plt.show()
