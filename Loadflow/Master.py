from Loadflow.PowerEq import *
import numpy as np
"""teta = np.zeros(3) #creates an array with the initial angles
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
"""
teta = np.array([0.0, 0.0, 0.0]) #creates an array with the initial angles
v = np.array([1.0, 1.0, 1.0]) #creates an array with the initial voltages
g12=(1/complex(0.05,0.2)).real
g13=(1/complex(0.05,0.1)).real
g23=(1/complex(0.05,0.15)).real
b12=(1/complex(0.05,0.2)).imag
b13=(1/complex(0.05,0.1)).imag
b23=(1/complex(0.05,0.15)).imag


g=np.array([[g12+g13,g12,g13],[g12,g12+g23,g23],[g13,g23,g13+g23]]) #creates a matrix with the conductanses
b=np.array([[b12+b13,b12,b13],[b12,b12+b23,b23],[b13,b23,b13+b23]]) #creates a matrix with the susceptanses
Pindex=np.array([1,2])#list with bus-numbers of the buses with known P
Qindex=np.array([1,2])#list with bus-numbers of the buses with known Q
vindex=np.array([1,2])#list with bus-numbers of the buses with unknown v
tindex=np.array([1,2])#list with bus-numbers of the buses with unknown angel
Pactual=np.array([-1.0,-0.5])
Qactual=np.array([-0.5,-0.5])
alpha = np.array([0.1,0.2])
beta = np.array([1,2])


#newtonrhapson_print(Pactual,Qactual,Pindex,Qindex,tindex,vindex,v,teta,g,b)
voltageStabilityFlatStart(Pactual,Qactual,Pindex, Qindex, tindex, vindex, g, b)
#voltageStabilityFlatStart(Pactual,Qactual,Pindex, Qindex, tindex, vindex, g, b)


