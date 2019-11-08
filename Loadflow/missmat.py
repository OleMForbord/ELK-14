import numpy as np
from Loadflow.powerinj import *

def missmat_p(Pactual,Pindex,v,teta,g,b):
    estpower=np.empty(0)
    for pi in range(0,Pindex.size):
        estpower=np.append(estpower,act_pow(Pindex[pi],v.size,v,teta,g,b))
    return Pactual-estpower

def missmat_q(Qactual,Qindex,v,teta,g,b):
    estpower=np.empty(0)
    for qi in range(0,Qindex.size):
        estpower=np.append(estpower,react_pow(Qindex[qi],v.size,v,teta,g,b))
    return Qactual-estpower

def power_missmatch(Pactual,Qactual,Pindex,Qindex,v,teta,g,b):
    missmat=np.append(missmat_p(Pactual,Pindex,v,teta,g,b),missmat_q(Qactual,Qindex,v,teta,g,b))

    return missmat