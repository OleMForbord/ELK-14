import sys
sys.path.append(".")
from Decoupled.fdlf import *
import numpy as np

def line_flow(i,j,teta,z, outage = False):
    if (outage == False):
        if z[i-1][j-1].imag!=0:
            return(teta[i-1] - teta[j-1]) / z[i-1][j-1].imag
        else:
            print('This is not a valid line')
            return 0
    else:
        return 0


def print_systemflow(teta, z,outage):
    for i in range(0, teta.size):
        for j in range(0, teta.size):
            if z[i][j].imag != 0:
                if outage:
                    print('flow between bus', i+1, j+1, ': 0')
                else:
                    print('flow between bus',i+1,j+1,': ',line_flow(i+1,j+1,teta,z))

def dcflow(Pactual,Pindex,tindex, v, teta, g, b,z):
    print('initial angles:\n',teta)
    b_matrix=heq_matrix(Pindex,tindex,z)
    print('B-matrix:\n',b_matrix)
    bmm_inv=np.linalg.inv(b_matrix)
    correction=bmm_inv.dot(missmat_p(Pactual,Pindex,v,teta,g,b))
    print('correction: ',correction)
    for a in range(0, tindex.size):
        teta[tindex[a] - 1] += correction[a]
    print('Updated angles: ',teta)
    pflow = np.zeros(shape=(teta.size, teta.size))
    for i in range(0, teta.size):
        for j in range(0, teta.size):
            if z[i][j].imag != 0:
                pflow[i][j] = line_flow(i+1,j+1,teta,z)
    print('power flow matrix:')
    return pflow
