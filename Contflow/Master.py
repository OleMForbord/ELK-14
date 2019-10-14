import numpy as np
from Contflow.ContinousPF import *

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
alpha = np.array([0,0])
beta = np.array([0.5,0.5])


contflow_print(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta,1)
