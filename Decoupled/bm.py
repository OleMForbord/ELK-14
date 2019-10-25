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

