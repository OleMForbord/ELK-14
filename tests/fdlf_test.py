from Loadflow.PowerEq import *
from Decoupled.FDLF import *
from Decoupled.optionally import *
import numpy as np

teta = np.array([0.0, 0.0, 0.0]) #creates an array with the initial angles
v = np.array([1.0, 1.0, 1.0]) #creates an array with the initial voltages
z12=(complex(0.05,0.2))
z13=complex(0.05,0.1)
z23=complex(0.05,0.15)

z=np.array([[0,z12,z13],[z12,0,z23],[z13,z23,0]])
g=g_matrix(z) #creates a matrix with the conductanses
b=b_matrix(z) #creates a matrix with the susceptances
x=z.imag

Pindex=np.array([1,2])#list with bus-numbers of the buses with known P
Qindex=np.array([1,2])#list with bus-numbers of the buses with known Q
vindex=np.array([1,2])#list with bus-numbers of the buses with unknown v
tindex=np.array([1,2])#list with bus-numbers of the buses with unknown angel
Pactual=np.array([-1.0,-0.5])
Qactual=np.array([-0.5,-0.5])
alpha = np.array([0.1,0.2])
beta = np.array([1,2])
#print_fdlf(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b,z,)
#print_standard(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b)
print(heq_matrix(Pindex,Qindex,tindex,vindex,v,teta,g,b))


