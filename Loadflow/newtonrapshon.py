import sys
sys.path.append(".")
from Loadflow.jacobi import *
import numpy as np
import matplotlib.pyplot as plt

def power_missmatch(Pactual,Qactual,Pindex,Qindex,v,teta,g,b):
    actualpower=np.append(Pactual,Qactual)

    estpower=np.empty(0)
    for pi in range(0,Pindex.size):
        estpower=np.append(estpower,act_pow(Pindex[pi],v.size,v,teta,g,b))
    for qi in range(0,Qindex.size):
        estpower=np.append(estpower,react_pow(Qindex[qi],v.size,v,teta,g,b))

    return actualpower-estpower

def newtonrhapson(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b):

    it = 0
    epsilonError = 0.001
    missmatch = power_missmatch(Pactual, Qactual, Pindex, Qindex, v, teta, g, b)
    jacobiinv = np.linalg.inv(jacobi(Pindex, Qindex, tindex, vindex, v, teta, g, b))
    correction = jacobiinv.dot(missmatch)
    print('correction: ',correction)
    print('teta: ',teta, 'correction[tindex.size: ',correction[:tindex.size])

    #from Master import CONVERGENCE_LIMIT

    while np.any(abs(correction) > epsilonError):
        if(it >= 1000):#CONVERGENCE_LIMIT):
            print("diverg")
            return 0
        for a in range(0, tindex.size):
            teta[tindex[a]-1] += correction[a]
        for vm in range(0, vindex.size):
            v[vindex[vm]-1] += correction[tindex.size + vm]

        missmatch = power_missmatch(Pactual, Qactual, Pindex, Qindex, v, teta, g, b)
        jacobiinv = np.linalg.inv(jacobi(Pindex, Qindex, tindex, vindex, v, teta, g, b))
        correction = jacobiinv.dot(missmatch)
        it += 1
    return 1

def newtonrhapson_print(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b):

    print('iteration number: 1')
    print('\nVoltage angles:', teta)
    print('Voltage magnitudes: ', v)
    it = 0
    epsilonError = 0.0001
    missmatch = power_missmatch(Pactual, Qactual, Pindex, Qindex, v, teta, g, b)
    print('Active power injection: ',pinj(v,teta,g,b))
    print('Reactive power injection: ', qinj(v, teta, g, b))
    jacobiinv = np.linalg.inv(jacobi(Pindex, Qindex, tindex, vindex, v, teta, g, b))
    correction = jacobiinv.dot(missmatch)

    print('\nJacobi-matrix:\n', jacobi(Pindex, Qindex, tindex, vindex, v, teta, g, b))
    print('\ncorrection: ', correction)
    #from Master import CONVERGENCE_LIMIT

    while np.any(abs(correction) > epsilonError)==True:
        if(it >= 1000):#CONVERGENCE_LIMIT):
            print("diverg")
            return 0
        for a in range(0, tindex.size):
            teta[tindex[a]-1] += correction[a]
        for vm in range(0, vindex.size):
            v[vindex[vm]-1] += correction[tindex.size + vm]
        it += 1

        missmatch = power_missmatch(Pactual, Qactual, Pindex, Qindex, v, teta, g, b)
        jacobiinv = np.linalg.inv(jacobi(Pindex, Qindex, tindex, vindex, v, teta, g, b))
        correction = jacobiinv.dot(missmatch)

        print('\niteration: ', it)
        print('\nVoltage angles:', teta)
        print('Voltage magnitudes: ', v)
        print('Active power injection: ', pinj(v, teta, g, b))
        print('Reactive power injection: ', qinj(v, teta, g, b))
        print('\nJacobi-matrix:\n', jacobi(Pindex, Qindex, tindex, vindex, v, teta, g, b))
        print('\ncorrection: ', correction)
    print("\nNumber of iterations before convergence is: ", it)
    for busi in vindex:
        print('Bus nr', busi, ': Voltage magnitude =', v[busi-1])


    plt.plot(-sum(Pactual), v[0], 'g^-', -sum(Pactual), v[1], 'bo-')


    return 1