import sys
sys.path.append(".")
from Decoupled.fdlf import *

def line_flow(i,j,teta,z):
    if z[i-1][j-1].imag!=0:
        return(teta[i-1] - teta[j-1]) / z[i-1][j-1].imag
    else:
        print('This is not a valid line')
        return 0



def dcflow(Pactual,Qactual,Pindex, Qindex, tindex, v, teta, g, b,z):
    print('initial angles:\n',teta)
    b_matrix=heq_matrix(Pindex,tindex,z)
    print('B-matrix:\n',b_matrix)
    bmm_inv=np.linalg.inv(b_matrix)
    correction=bmm_inv.dot(missmat_p(Pactual,Qactual,Pindex,Qindex,tindex,v,teta,g,b))
    print('correction: ',correction)
    for a in range(0, tindex.size):
        teta[tindex[a] - 1] += correction[a]
    print('Updated angles: ',teta)
    pflow = np.zeros(shape=(teta.size, teta.size))
    for i in range(0, teta.size):
        for j in range(0, teta.size):
            if z[i][j] != 0:
                pflow[i][j] = line_flow(i+1,j+1,teta,z)
    return pflow
