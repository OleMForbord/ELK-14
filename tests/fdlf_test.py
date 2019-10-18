import sys
sys.path.append(".")
from tests.basecase import *
from Decoupled.fdlf import *

#print_dual(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b,z)
#print_primal(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b,z)
print_standard(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b,z)



