import sys
sys.path.append(".")
from Decoupled.bmm import *


def distfactor(Pindex,tindex,z,i,j,slacknr):
    if i==slacknr or j==slacknr:
        #if the line is connectet to the slack bus, all ptdf is 1, except for slack bus which is 0
        ptdf=np.ones(tindex.size+1)
        ptdf[slacknr-1]=0
        return ptdf
    heq=heq_matrix(Pindex,tindex,z)
    x_vec = np.zeros(heq[0].size)
    x_vec[i-1] = 1/(z[i-1][j-1].imag)
    x_vec[j-1] = -x_vec[i-1]
    ptdf=np.linalg.inv(heq).dot(x_vec)
    ptdf=np.insert(ptdf,slacknr-1,0)
    return ptdf

