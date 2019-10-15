from Loadflow.PowerEq import *
from Decoupled.FDLF import *
from Decoupled.optionally import *
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
z12=(complex(0.05,0.2))
z13=complex(0.05,0.1)
z23=complex(0.05,0.15)

z=np.array([[0,z12,z13],[z12,0,z23],[z13,z23,0]])
g=g_matrix(z) #creates a matrix with the conductanses
b=b_matrix(z) #creates a matrix with the susceptances
x=z.imag
print(x)
Pindex=np.array([1,2])#list with bus-numbers of the buses with known P
Qindex=np.array([1,2])#list with bus-numbers of the buses with known Q
vindex=np.array([1,2])#list with bus-numbers of the buses with unknown v
tindex=np.array([1,2])#list with bus-numbers of the buses with unknown angel
Pactual=np.array([-1.0,-0.5])
Qactual=np.array([-0.5,-0.5])
alpha = np.array([0.1,0.2])
beta = np.array([1,2])

print()
#print_fdlf(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b,z,)
#print_primal(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b,z)
#print(h_matrix(Pindex,tindex,v,teta,g,b))
#print(b)
#newtonrhapson_print(Pactual,Qactual,Pindex,Qindex,tindex,vindex,v,teta,g,b)
#voltageStabilityFlatStart(Pactual,Qactual,Pindex, Qindex, tindex, vindex, g, b)
#voltageStabilityFlatStart(Pactual,Qactual,Pindex, Qindex, tindex, vindex, g, b)


