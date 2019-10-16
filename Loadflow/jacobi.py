import sys
sys.path.append(".")
from Loadflow.powerinj import *

def jacobi(Pindex, Qindex, tindex, vindex, v, teta, g, b):
    jac = np.empty(shape=(Pindex.size+ Qindex.size, tindex.size + vindex.size))
    for pi in range(0,Pindex.size):
        row = np.empty(0)
        row = np.append(row, dpdt(Pindex[pi], tindex, v, teta, g, b))
        row = np.append(row, dpdv(Pindex[pi], vindex, v, teta, g, b))
        jac[pi] = row
    for qi in range(0,Qindex.size):
        row = np.empty(0)
        row = np.append(row, dqdt(Qindex[qi], tindex, v, teta, g, b))
        row = np.append(row, dqdv(Qindex[qi], vindex, v, teta, g, b))
        jac[Pindex.size + qi] = row
    return jac