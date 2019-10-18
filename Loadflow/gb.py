import sys
sys.path.append(".")
import numpy as np

def g_matrix(z):
    g=np.zeros(shape=(z[0].size,z[0].size))
    for i in range (0,z[0].size):
        for j in range (0,z[0].size):
            if j!=i:
                g[i][j]=(1/z[i][j]).real
        g[i][i]=np.sum(g[i])
    return g
def b_matrix(z):
    b=np.zeros(shape=(z[0].size,z[0].size))
    for i in range (0,z[0].size):
        for j in range (0,z[0].size):
            if j!=i:
                b[i][j]=(1/z[i][j]).imag
        b[i][i]=np.sum(b[i])
    return b