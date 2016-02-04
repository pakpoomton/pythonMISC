import numpy as np
import matplotlib.pyplot as plt

N = 10000 #number of term in the series
Z = 100 # number of point in space to consider

MM = 10000
A = MM*1.
B = MM*0.01
D = 1
w = 1.
L = 4.
pi = np.pi
t = 1000 # time till the end


x = np.linspace(0,L,100)
U = np.ones(Z)
U = U + (1-np.exp(-B*t))*A*w/(L*B)



for n in range(1,N+1):
   coeff = (D*n*n*pi*pi/L/L) + B
   U = U + (2*A/(n*pi*coeff))*(1-np.exp(-coeff*t))*np.sin(n*pi*w/L)*np.cos(n*pi*x/L)
   print(str(n))
   
#print(str(U))
plt.plot(x, U, color = 'k')
plt.xlabel('location')
plt.ylabel('temperature')
plt.axis([0, L, -0.1, 1.1*100])
plt.show()




