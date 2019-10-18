import sys
sys.path.append(".")
import numpy as np

def h_matrix(Pindex,tindex,b):
    matrix= np.zeros(shape=(Pindex.size,tindex.size))
    for i in Pindex:
        row=np.zeros(tindex.size)
        for j in tindex:
            if(j!=i):
                row[j-1]=b[i-1][j-1]
            else:
                row[j - 1] = -b[i - 1][j - 1]
        matrix[i-1]=row
    return matrix

def heq_matrix(Pindex,tindex,z):
    matrix= np.zeros(shape=(Pindex.size,tindex.size))
    for i in Pindex:
        for j in tindex:
            if j!=i:
                matrix[i-1][j-1] = -1/(z[i-1][j-1].imag)
        for k in range(0, z[i-1].size):
            if z[i-1][k] != 0:
                matrix[i-1][i-1] += np.sum(1 / (z[i-1][k].imag))
    return matrix

