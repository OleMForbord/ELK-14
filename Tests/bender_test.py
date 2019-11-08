import sys
sys.path.append(".")
from Loadflow.voltstab import *
from DCflow.dcflow import *
from Contingency.imml import *
import numpy as np
from Contingency.econ_disp import *

teta = np.zeros(4) #creates an array with the initial angles
v = np.array([1.0, 1.0, 1.0]) #creates an array with the initial voltages
z12=(complex(0.0,0.2))
z13=complex(0.0,0.1)
z23=complex(0.0,0.25)
z34=complex(0.0,0.25)

z=np.array([[0,z12,z13,inf],[z12,0,z23,inf],[z13,z23,0,z34],[inf,inf,z34,0]])
g=g_matrix(z) #creates a matrix with the conductanses
b=b_matrix(z) #creates a matrix with the susceptances
x=z.imag
econ_disp()