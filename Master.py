import np as np

from PowerEq import *
import numpy as np

teta = np.zeros(3) #creates an array with the initial angles
#v = np.array([1.0, 1.0, 1.0]) #creates an array with the initial voltages
v = np.full((1,3), 1, dtype = float)

r = np.array([0.1, 0.05, 0.05])
x = np.array([0.2, 0.1, 0.15])

z = r + 1j*x
y = 1/z
g = np.array([[y[0].real+y[1].real, y[0].real, y[1].real], [y[0].real, y[0].real + y[1].real, y[2].real], [y[1].real, y[2].real, y[1].real + y[2].real]])
j = 0
g = np.zeros((3,3))
for i in range(len(y)):
    if(i == j):
        g[i][j] = y[i].real+y[j+1].real



#g=np.array([[g12+g13,g12,g13],[g12,g12+g23,g23],[g13,g23,g13+g23]]) #creates a matrix with the conductanses
#b=np.array([[b12+b13,b12,b13],[b12,b12+b23,b23],[b13,b23,b13+b23]]) #creates a matrix with the susceptanses
Pindex=np.array([1,2])#list with bus number of buses with known P
Qindex=np.array([1,2])#list with bus number of buses with known Q
vindex=np.array([1,2])#list with bus number of buses with unknown v
tindex=np.array([1,2])#list with bus number of buses with unknown angel
Pactual=np.array([-1.0,-0.5])
Qactual=np.array([-0.5,-0.5])


newtonrhapson(Pactual,Qactual,Pindex,Qindex,tindex,vindex,v,teta,g,b,0.01)





