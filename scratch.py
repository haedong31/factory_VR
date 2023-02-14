import numpy as np
import matplotlib.pyplot as plt

def learning_curve(a,b,x):
    return a*np.power(x,-b)

r1 = 0.8
r2 = 0.6
a = 1
b1 = -(np.log(r1)/np.log(2))
b2 = -(np.log(r2)/np.log(2))

x = np.arange(1,100)

y1 = learning_curve(a, b1, x)
y2 = learning_curve(a, b2, x)

plt.plot(x,y1)
plt.plot(x,y2)
plt.xlabel('Number of units produced')
plt.ylabel('Production time')
plt.legend(['80% Learning curve', '60% Learning curve'])