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
    #print(actualpower)
    estpower=np.empty(0)
    for pi in range(0,Pindex.size):
        estpower=np.append(estpower,act_pow(Pindex[pi],v.size,v,teta,g,b))
    for qi in range(0,Qindex.size):
        estpower=np.append(estpower,react_pow(Qindex[qi],v.size,v,teta,g,b))
    #print(estpower)
    return actualpower-estpower

def newtonrhapson(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b):
    it = 1
    epsilonError = 0.001
    deltapower = power_missmatch(Pactual, Qactual, Pindex, Qindex, v, teta, g, b)
    jacobiinv = np.linalg.inv(jacobi(Pindex, Qindex, tindex, vindex, v, teta, g, b))
    correction = jacobiinv.dot(deltapower)
    from Master import CONVERGENCE_LIMIT
    #while abs(correction[2]) > epsilonError:
    while abs(correction[2]) > epsilonError:
        if(it >= CONVERGENCE_LIMIT):

            return 0

        #print('Jacobi-matrix:\n', jacobi(Pindex, Qindex, tindex, vindex, v, teta, g, b))
        i = 0
        #print("\nCorection: \n", correction)
        for a in range(0, tindex.size):
            teta[a] = teta[a] + correction[a]
            #print('Bus nr', tindex[a], ': Voltage angle =', teta[a])

        for vm in range(0, vindex.size):

            v[vm] = v[vm] + correction[tindex.size - 1 + vm]
            #print('Bus nr', vindex[vm], ': Voltage magnitude =', v[vm])
        deltapower = power_missmatch(Pactual, Qactual, Pindex, Qindex, v, teta, g, b)
        jacobiinv = np.linalg.inv(jacobi(Pindex, Qindex, tindex, vindex, v, teta, g, b))
        correction = jacobiinv.dot(deltapower)
        it += 1
    print("Number of iterations before convergence is: ", it)
    print('Bus nr', vindex[0], ': Voltage magnitude =', v[0])
    print('Bus nr ', vindex[1], ': Voltage magnitude =', v[1])

    plt.plot(-sum(Pactual), v[0], 'g^-', -sum(Pactual), v[1], 'bo-')


    return 1
    #for i in correction:
       # while i>0.5:
        #    newtonrhapson(Pactual,Qactual,Pindex,Qindex,tindex,vindex,v,teta,g,b)
