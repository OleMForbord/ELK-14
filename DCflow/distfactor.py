import sys
sys.path.append(".")
from Decoupled.bm import *
<<<<<<< HEAD
=======

>>>>>>> 437774f9083faeb3c232c6b8b5ae35f465196aa8

def distfactor(Pindex,tindex,z,i,j):

    heq=heq_matrix(Pindex,tindex,z)
    x_vec = np.zeros(heq[0].size)
    x_vec[i-1] = 1/(z[i-1][j-1].imag)
    x_vec[j-1] = -x_vec[i-1]

    return np.linalg.inv(heq).dot(x_vec)

