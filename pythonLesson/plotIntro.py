import matplotlib.pyplot as plt

time = [1,2,3,4,5,6,7,8,9,10]
population1 = [10,20,30,40,50,60,70,80,90,100]
population2 = [1,2,4,8,16,32,64,128,256,512]


plt.plot(time, population1, color = 'r', label='linear')
plt.plot(time, population2, color = 'b', label = 'exponent')

plt.title('Population Growth')
plt.legend(loc = 2) # location of figure legend
plt.xlabel('Time')
plt.ylabel('Population')
plt.show()
