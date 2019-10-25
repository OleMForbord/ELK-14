import sys
sys.path.append(".")
from Decoupled.bm import *

def distfactor(Pindex,tindex,z,i,j):

    heq=heq_matrix(Pindex,tindex,z)
    heq_new = np.array([heq[i-1],heq[j-1]])
    print(heq_new)
    return np.linalg.inv(heq_new).dot([1/z[i-1][j-1].imag,-1/z[i-1][j-1].imag])