import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

# set parameters of PDE
a = 2.8e-4
b = 5e-3
tau = .1
k = -.005

size = 100 # size of the 2D grid
# discretized to a matrix of size x size
dx = 2./size  # space step    ..domain size = 2io
# set parameters for discretize time and space
# note that: dt<=dx**2/2 ensures stable scheme

timeTotal = 10.0  # total time
dt = .9 * dx**2/2  # time step
n = int(timeTotal/dt)

# u  represent concentration of a substance favoring
# skin pigmantation; v represents another substance
# interacting with u and prevent pigmentation. Here
# we start from uniform random u and v each grid point in space
U = np.random.rand(size, size)
V = np.random.rand(size, size)

def laplacian(Z): #discrete estimation of dZ matrix
    Ztop = Z[0:-2,1:-1]
    Zleft = Z[1:-1,0:-2]
    Zbottom = Z[2:,1:-1]
    Zright = Z[1:-1,2:]
    Zcenter = Z[1:-1,1:-1]
    return (Ztop + Zleft + Zbottom + Zright - 4 * Zcenter) / dx**2




outputMatU = []
outputMatU.append(sp.copy(U))

# We simulate the PDE with the finite difference method.
for i in range(n):
# We compute the Laplacian of u and v...used as estimated changes
# of U and V in each time step
    deltaU = laplacian(U)
    deltaV = laplacian(V)
    # We take the values of u and v inside the grid.
    Uc = U[1:-1,1:-1]
    Vc = V[1:-1,1:-1]
    # We update the variables.
    U[1:-1,1:-1], V[1:-1,1:-1] = \
        Uc + dt * (a * deltaU + Uc - Uc**3 - Vc + k), \
        Vc + dt * (b * deltaV + Uc - Vc) / tau
    # Neumann conditions: derivatives at the edges
    # are null.
    for Z in (U, V):
        Z[0,:] = Z[1,:]
        Z[-1,:] = Z[-2,:]
        Z[:,0] = Z[:,1]
        Z[:,-1] = Z[:,-2]

    outputMatU.append(sp.copy(U))

#Uf = outputMatU[len(outputMatU)-1]
    
#plt.imshow(Uf, cmap=plt.cm.hot, extent=[-1,1,-1,1]);
#plt.xticks([]); plt.yticks([]);
#plt.show()

#display result. We have 6 subfigure showing results from different diffusion time points. 
intv_ind = int(len(outputMatU)/5) # interval index between each plot
intv_time = timeTotal/5 # interval time between each plot

f, axarr = plt.subplots(2, 3)
u = outputMatU[0]
axarr[0, 0].imshow( u, cmap=plt.cm.hot, interpolation='nearest', origin='lower')
axarr[0, 0].set_title('time = 0')
u = outputMatU[intv_ind - 1]
axarr[0, 1].imshow(u, cmap=plt.cm.hot, interpolation='nearest', origin='lower')
axarr[0, 1].set_title('time = ' + str(intv_time))
u = outputMatU[2*intv_ind - 1]
axarr[0, 2].imshow(u, cmap=plt.cm.hot, interpolation='nearest', origin='lower')
axarr[0, 2].set_title('time = ' + str(2*intv_time))
u = outputMatU[3*intv_ind - 1]
axarr[1, 0].imshow(u, cmap=plt.cm.hot, interpolation='nearest', origin='lower')
axarr[1, 0].set_title('time = ' + str(3*intv_time))
u = outputMatU[4*intv_ind - 1]
axarr[1, 1].imshow(u, cmap=plt.cm.hot, interpolation='nearest', origin='lower')
axarr[1, 1].set_title('time = ' + str(4*intv_time))
u = outputMatU[5*intv_ind - 1]
axarr[1, 2].imshow(u, cmap=plt.cm.hot, interpolation='nearest', origin='lower')
axarr[1, 2].set_title('time = ' + str(5*intv_time))

plt.show()
    




