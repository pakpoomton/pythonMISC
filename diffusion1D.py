from fipy import *

# define 1D domain in which we will divide into a grid for our numeric solver
nx = 50 # number of solution points
dx = 1. # grid spacing
mesh = Grid1D(nx=nx, dx=dx) # define an obj represent linear structure grid
x = mesh.cellCenters[0] # ??
L = nx * dx # domain length

phi = CellVariable(mesh=mesh, name=r"$\phi$")

phi.value = 0.
phi.setValue(1., where=(x > L/2. - L/10.) & (x < L/2. + L/10.))
viewer = Viewer(vars=phi, datamin=-0.1, datamax=1.1)

D = 1.

eq = TransientTerm() == DiffusionTerm(D)
dt = 10. * dx**2 / (2 * D)

steps = 200
for step in range(steps):
       eq.solve(var=phi, dt=dt)
       viewer.plot()


