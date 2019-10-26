import sys
sys.path.append(".")
from Decoupled.bmm import *

def angles_imml(i,j,delta_h,Pactual,Pindex,tindex,z):
    i=np.where(Pindex==i)
    j=np.where(Pindex==j)
    if i[0].size==0 or j[0].size==0:
        print('One of the bus you entered is the slack bus, invalid operation.')
        return 0
    i=i[0][0]
    j=j[0][0]
    h_inv=np.linalg.inv(heq_matrix(Pindex,tindex,z))
    teta_init=h_inv.dot(Pactual)
    c=1/(1/delta_h+h_inv[i][i]+h_inv[j][j]-h_inv[i][j]-h_inv[j][i])
    return teta_init-(h_inv[:,i]-h_inv[:,j])*c*(teta_init[i]-teta_init[j])
