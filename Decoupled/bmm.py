import sys
sys.path.append(".")
import numpy as np

def l_matrix(Qindex,vindex,b):
    matrix= np.zeros(shape=(Qindex.size,vindex.size))
    for i in Qindex:
        row=np.zeros(vindex.size)
        for j in vindex:
            if(j!=i):
                row[j-1]=b[i-1][j-1]
            else:
                row[j - 1] = -b[i - 1][j - 1]
        matrix[i-1]=row
    return matrix
    
def leq_matrix(Qindex,vindex,z):
    matrix= np.zeros(shape=(Qindex.size,vindex.size))
    for i in Qindex:
        for j in vindex:
            if j!=i:
                matrix[i-1][j-1] = -1/(z[i-1][j-1].imag)
        for k in range(0, z[i-1].size):
            if z[i-1][k] != 0:
                matrix[i-1][i-1] += np.sum(1 / (z[i-1][k].imag))
    return matrix
