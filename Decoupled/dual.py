import sys
sys.path.append(".")
from Decoupled.bm import *
from Decoupled.bmm import *
from Loadflow.missmat import *

def print_dual(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b,z):
    it = 1
    epsilonError = 0.001
    l_inv=np.linalg.inv(l_matrix(Qindex,vindex,b))
    heq_inv=np.linalg.inv(heq_matrix(Pindex,tindex,z))

    correction_teta = heq_inv.dot(missmat_p(Pactual, Pindex, v, teta, g, b))
    correction_v = l_inv.dot(missmat_q(Qactual, Qindex, v, teta, g, b))


    print('l-matrix:\n',np.linalg.inv(l_inv))
    print('\nheq_marix:\n',np.linalg.inv(heq_inv))
    print('\nVoltage angles: ', teta)
    print('Voltage magnetudes: ', v)

    while np.any(abs(correction_teta)  > epsilonError):
        if it>1000:
            print('No solution found')
            return 0
        correction_v = l_inv.dot(missmat_q( Qactual,Qindex, v, teta, g, b))
        for vm in range(0, vindex.size):
            v[vindex[vm] - 1] += correction_v[vm]
        print('\nMagnitudes correction:', correction_v)
        print('Voltage magnetudes: ', v)

        correction_teta = heq_inv.dot(missmat_p(Pactual, Pindex, v, teta, g, b))
        for a in range(0,tindex.size):
            teta[tindex[a]-1] += correction_teta[a]
        print('\niteraton: ', it)
        print('\nAngle correction:', correction_teta)
        print('Voltage angles: ', teta)
        it+=1
    return 1
