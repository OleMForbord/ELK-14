import sys
sys.path.append(".")
from Decoupled.bm import *
from Decoupled.bmm import *
from Decoupled.mismat import *

def print_primal(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b,z):
    it = 1
    epsilonError = 0.001
    h_inv=np.linalg.inv(h_matrix(Pindex,tindex,b))
    leq_inv=np.linalg.inv(leq_matrix(Qindex,vindex,z))

    correction_teta = h_inv.dot(missmat_p(Pactual, Qactual, Pindex, Qindex, tindex, v, teta, g, b))
    correction_v = leq_inv.dot(missmat_q(Pactual, Qactual, Pindex, Qindex, vindex, v, teta, g, b))


    print('h-matrix:\n',np.linalg.inv(h_inv))
    print('\nleq-matrix:\n',np.linalg.inv(leq_inv))
    print('\nVoltage angles: ', teta)
    print('Voltage magnetudes: ', v)

    while np.any(abs(correction_teta)  > epsilonError):
        if(it>1000):
            print('the primal method did not converge, trying the dual:')
            return 0
        correction_teta = h_inv.dot(missmat_p(Pactual, Qactual, Pindex, Qindex, tindex, v, teta, g, b))
        for a in range(0,tindex.size):
            teta[tindex[a]-1] += correction_teta[a]
        print('\niteraton: ', it)
        print('\nAngle correction:', correction_teta)
        print('Voltage angles: ', teta)

        correction_v = leq_inv.dot(missmat_q(Pactual, Qactual, Pindex, Qindex, vindex, v, teta, g, b))
        for vm in range(0, vindex.size):
            v[vindex[vm]-1] += correction_v[vm]
        print('\nMagnitudes correction:',correction_v)
        print('Voltage magnetudes: ', v)

        it += 1
    it = 1
    return 1
