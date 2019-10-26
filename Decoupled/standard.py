import sys
sys.path.append(".")
from Decoupled.bmm import *
from Loadflow.missmat import *

def print_standard(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b,z):
    it = 1
    epsilonError = 0.001
    matrix=np.zeros(shape=(Pindex.size+Qindex.size,tindex.size+vindex.size))
    matrix[:Pindex.size,0:tindex.size]=heq_matrix(Pindex,tindex,z)
    matrix[Pindex.size:,tindex.size:] = leq_matrix(Qindex,vindex,z)

    matrix_inv=np.linalg.inv(matrix)
    correction=matrix_inv.dot(power_missmatch(Pactual,Qactual,Pindex,Qindex,v,teta,g,b))
    print('standard decoupled matrix:\n',matrix)
    print('\nVoltage angles: ', teta)
    print('Voltage magnetudes: ', v)

    while np.any(abs(correction)  > epsilonError):
        if(it>1000):
            print('the primal method did not converge, trying the dual:')
            return 0
        correction = matrix_inv.dot(power_missmatch(Pactual, Qactual, Pindex, Qindex, v, teta, g,b))
        for a in range(0,tindex.size):
            teta[tindex[a]-1] += correction[a]
        for vm in range(0, vindex.size):
            v[vindex[vm]-1] += correction[tindex.size+vm]
        print('\niteraton: ', it)
        print('\nAngle correction:', correction[Pindex.size:])
        print('Voltage angles: ', teta)
        print('\nMagnitudes correction:',correction[:Pindex.size])
        print('Voltage magnetudes: ', v)




        correction= matrix_inv.dot(power_missmatch(Pactual, Qactual, Pindex, Qindex, v, teta, g, b))

        it += 1
    it = 1
    return 1
