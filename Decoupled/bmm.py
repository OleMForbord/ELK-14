import sys
sys.path.append(".")
import numpy as np

def heq_matrix(Pindex,tindex,z):
    matrix= np.zeros(shape=(Pindex.size,tindex.size))
    for i in Pindex:
        for j in tindex:
            if j!=i:
                matrix[i-1][j-1] = -1/(z[i-1][j-1].imag)
        for k in range(0, z[i-1].size):
            if z[i-1][k].imag != 0:
                matrix[i-1][i-1] += np.sum(1 / (z[i-1][k].imag))
    return matrix
    
def leq_matrix(Qindex,vindex,z):
    matrix= np.zeros(shape=(Qindex.size,vindex.size))
    for i in Qindex:
        for j in vindex:
            if j!=i:
                matrix[i-1][j-1] = -1/(z[i-1][j-1].imag)
        for k in range(0, z[i-1].size):
            if z[i-1][k].imag != 0:
                matrix[i-1][i-1] += np.sum(1 / (z[i-1][k].imag))
    return matrix
