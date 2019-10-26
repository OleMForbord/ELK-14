import sys
sys.path.append(".")
from Decoupled.primal import *
from Decoupled.dual import *
from Decoupled.standard import *

def fdlf_print(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b,z):
    if(print_primal(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b,z)==0):
        print_dual(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b,z)

