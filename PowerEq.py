from math import *
import numpy as np
import matplotlib.pyplot as plt

def tij(
    i, j, teta, g, b
):  # Creates the tij needed in the derivatives. i,j=bus nr i and j.teta=the angles g=conductance matrix, b=susceptance matrix
    return g[i - 1][j - 1] * cos(teta[i - 1] - teta[j - 1]) + b[i - 1][j - 1] * sin(
        teta[i - 1] - teta[j - 1]
    )


def uij(i, j, teta, g, b):  # Creeates the uij needed in the derivatives
    return g[i - 1][j - 1] * sin(teta[i - 1] - teta[j - 1]) - b[i - 1][j - 1] * cos(
        teta[i - 1] - teta[j - 1]
    )


def act_pow(busi, n, v, teta, g, b):  # n=number of buses
    p = (v[busi-1]**2) * g[busi-1][busi-1]
    for busj in range(1, n+1):
        if busj == busi:
            continue
        else:
            p = p - v[busi-1] * v[busj-1] * tij(busi, busj, teta, g, b)
    return p
    """pi = 0
    for busj in range(1, n+1):

        pi += (v[busi-1] * v[busj-1] * tij(busi, busj, teta, g,b))

    return pi"""



def react_pow(busi, n, v, teta, g, b):  # n=number of buses
   q = -(v[busi-1]**2) * b[busi-1][busi-1]
   for busj in range(1, n + 1):
       if busj == busi:
           continue
       else:
           q = q - v[busi-1] * v[busj-1] * uij(busi, busj, teta, g, b)
   return q

   """qi = 0
   for busj in range(1, n + 1):
        qi += (v[busi-1] * v[busj-1] * uij(busi, busj, teta, g, b))
   return qi"""

def dpdt(
    busi, tindex, v, teta, g, b
):
    row = np.empty(0)
    for busj in tindex:
        if busj == busi:
            element = 0
            for vj in range(1,v.size+1):
                if vj != busi:
                    element += v[busi-1] * v[vj-1] * uij(busj, vj, teta, g, b)
            row = np.append(row, element)
        else:
            row = np.append(
                row, -v[busi-1] * v[busj-1] * uij(busi, busj, teta, g, b)
            )
    return row


def dpdv(busi, vindex, v, teta, g, b):
    row = np.empty(0)
    for busj in vindex:
        if busj==busi:
            element = 0
            for vj in range(1,v.size+1):
                if vj==busj:
                    element += 2 * v[busi-1] * g[busi-1][busi-1]
                else:
                    element += -v[vj-1] * tij(busi, vj, teta, g, b)
            row = np.append(row, element)
        else:
            row = np.append(row, -v[busj-1] * tij(busi, busj, teta, g, b))
    return row


def dqdt(busi, tindex, vindex, v, teta, g, b):
    row = np.empty(0)
    for busj in tindex:
        if busj==busi:
            element = 0
            for vj in range(1,vindex.size+1):
                if vj != busi:
                    element += (
                        -v[busi-1] * v[vj-1] * tij(busi, busj, teta, g, b)
                    )
            row = np.append(row, element)
        else:
            for vj in range(1,vindex.size+1):
                if vj!=busi:
                    row = np.append(
                        row, v[busi-1] * v[vj-1] * tij(busi, busj, teta, g, b)
            )
    return row


def dqdv(busi, vindex, v, teta, g, b):
    row = np.empty(0)
    for busj in vindex:
        if busj == busi:
            element = 0
            for vj in range(1,vindex.size+1):
                if vj==busj:
                    element += -2 * v[busi-1] * b[busi-1][busi-1]
                else:
                    element += -v[vj-1] * uij(busi, busj, teta, g, b)
            row = np.append(row, element)
        else:
            row = np.append(row, -v[busi-1] * uij(busi, busj, teta, g, b))
    return row


def jacobi(Pindex, Qindex, tindex, vindex, v, teta, g, b):
    jac = np.empty(shape=(Pindex.size+ Qindex.size, tindex.size + vindex.size))
    for pi in range(0,Pindex.size):
        row = np.empty(0)
        row = np.append(row, dpdt(Pindex[pi], tindex, v, teta, g, b))
        row = np.append(row, dpdv(Pindex[pi], vindex, v, teta, g, b))
        jac[pi] = row
    for qi in range(0,Qindex.size):
        row = np.empty(0)
        row = np.append(row, dqdt(Qindex[qi], tindex, vindex, v, teta, g, b))
        row = np.append(row, dqdv(Qindex[qi], vindex, v, teta, g, b))
        jac[Pindex.size + qi] = row
    return jac

def power_missmatch(Pactual,Qactual,Pindex,Qindex,v,teta,g,b):
    actualpower=np.append(Pactual,Qactual)

    estpower=np.empty(0)
    for pi in range(0,Pindex.size):
        estpower=np.append(estpower,act_pow(Pindex[pi],v.size,v,teta,g,b))
    for qi in range(0,Qindex.size):
        estpower=np.append(estpower,react_pow(Qindex[qi],v.size,v,teta,g,b))

    return actualpower-estpower

def missmatch_print(Pactual,Qactual,Pindex,Qindex,v,teta,g,b):
    actualpower=np.append(Pactual,Qactual)

    estpower=np.empty(0)
    for pi in range(0,Pindex.size):
        estpower=np.append(estpower,act_pow(Pindex[pi],v.size,v,teta,g,b))
    print('\nActive power injection: ', estpower)
    for qi in range(0,Qindex.size):
        estpower=np.append(estpower,react_pow(Qindex[qi],v.size,v,teta,g,b))
    print('Reactive power injection: ',estpower[:Pindex.size])
    missmatch=actualpower-estpower

    print('\nActive power missmatch: ', missmatch[:Pindex.size])
    print('Reactive power missmatch: ', missmatch[Pindex.size:])
    return actualpower-estpower

def voltageStabilityFlatStart(Pactual,Qactual,Pindex, Qindex, tindex, vindex, g, b):
    teta = np.array([0.0, 0.0, 0.0])  # creates an array with the initial angles
    v = np.array([1.0, 1.0, 1.0])  # creates an array with the initial voltages
    syst_load = abs(Pactual[0] + Pactual[1])
    P_load = np.array([syst_load])

    newtonrhapson(Pactual, Qactual, Pindex, Qindex, tindex, vindex, v, teta, g, b)
    v1 = np.array(v[0])
    v2 = np.array(v[1])

    while True:
        teta = np.array([0.0, 0.0, 0.0])  # creates an array with the initial angles
        v = np.array([1.0, 1.0, 1.0])  # creates an array with the initial voltages
        if (newtonrhapson(Pactual, Qactual, Pindex, Qindex, tindex, vindex, v, teta, g, b) == 0):
            break
        Pactual = Pactual - [0.2 * 0.3, 0.2 * 0.7]
        # Qactual = Qactual - [0.2 * 0.3, 0.2 * 0.7]
        syst_load += abs(Pactual[0] + Pactual[1])
        newtonrhapson(Pactual, Qactual, Pindex, Qindex, tindex, vindex, v, teta, g, b)
        P_load = np.append(P_load, syst_load)
        v1 = np.append(v1, v[0])
        v2 = np.append(v2, v[1])
    print("Divergence at Pactual: ", Pactual)
    print(v1, v2, P_load)

    plt.xlabel("Load power")
    plt.ylabel("Voltage")
    plt.plot(P_load, v1)
    plt.plot(P_load, v2)
    plt.show()

def voltageStabilityAccumulating(Pactual,Qactual,Pindex, Qindex, tindex, vindex, teta, v,  g, b):
    while True:
        print("Teta is:", teta)
        if (newtonrhapson(Pactual, Qactual, Pindex, Qindex, tindex, vindex, v, teta, g, b) == 0):
            break
            Pactual = Pactual - [0.2 * 0.3, 0.2 * 0.7]
        # Qactual = Qactual - [0.2 * 0.3, 0.2 * 0.7]

    print("Divergence at Pactual: ", Pactual)
    plt.xlabel("Load power")
    plt.ylabel("Voltage")
    plt.show()




def newtonrhapson(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b):

    it = 1
    epsilonError = 0.001
    missmatch = power_missmatch(Pactual, Qactual, Pindex, Qindex, v, teta, g, b)
    jacobiinv = np.linalg.inv(jacobi(Pindex, Qindex, tindex, vindex, v, teta, g, b))
    correction = jacobiinv.dot(missmatch)

    #from Master import CONVERGENCE_LIMIT

    while abs(correction[2]) > epsilonError:
        if(it >= 1000):#CONVERGENCE_LIMIT):
            print("diverg")
            return 0
        for a in range(0, tindex.size):
            teta[a] = teta[a] + correction[a]
        for vm in range(0, vindex.size):
            v[vm] = v[vm] + correction[tindex.size + vm]

        missmatch = power_missmatch(Pactual, Qactual, Pindex, Qindex, v, teta, g, b)
        jacobiinv = np.linalg.inv(jacobi(Pindex, Qindex, tindex, vindex, v, teta, g, b))
        correction = jacobiinv.dot(missmatch)
        it += 1
    return 1

def newtonrhapson_print(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b):

    print('iteration number: 1')
    print('\nVoltage angles:', teta)
    print('Voltage magnitudes: ', v)
    it = 1
    epsilonError = 0.001
    missmatch = missmatch_print(Pactual, Qactual, Pindex, Qindex, v, teta, g, b)
    jacobiinv = np.linalg.inv(jacobi(Pindex, Qindex, tindex, vindex, v, teta, g, b))
    correction = jacobiinv.dot(missmatch)

    print('\nJacobi-matrix:\n', jacobi(Pindex, Qindex, tindex, vindex, v, teta, g, b))
    print('\ncorrection: ', correction)
    #from Master import CONVERGENCE_LIMIT

    while abs(correction[2]) > epsilonError:
        if(it >= 1000):#CONVERGENCE_LIMIT):
            print("diverg")
            return 0
        for a in range(0, tindex.size):
            teta[a] = teta[a] + correction[a]
        for vm in range(0, vindex.size):
            v[vm] = v[vm] + correction[tindex.size + vm]
        it += 1
        print('\niteration: ', it)
        print('\nVoltage angles:', teta)
        print('Voltage magnitudes: ', v)
        missmatch = missmatch_print(Pactual, Qactual, Pindex, Qindex, v, teta, g, b)
        jacobiinv = np.linalg.inv(jacobi(Pindex, Qindex, tindex, vindex, v, teta, g, b))
        correction = jacobiinv.dot(missmatch)

        print('\nJacobi-matrix:\n', jacobi(Pindex, Qindex, tindex, vindex, v, teta, g, b))
        print('\ncorrection: ', correction)
    print("\nNumber of iterations before convergence is: ", it)
    print('Bus nr', vindex[0], ': Voltage magnitude =', v[0])
    print('Bus nr ', vindex[1], ': Voltage magnitude =', v[1])

    plt.plot(-sum(Pactual), v[0], 'g^-', -sum(Pactual), v[1], 'bo-')


    return 1