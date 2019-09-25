from math import *
import numpy as np

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


def react_pow(busi, n, v, teta, g, b):  # n=number of buses
    q = -(v[busi-1]**2) * b[busi-1][busi-1]
    for busj in range(1, n+1):
        if busj == busi:
            continue
        else:
            print('it did really continue')
            q = q - v[busi-1] * v[busj-1] * uij(busi, busj, teta, g, b)
    return q


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
    jac = np.empty(shape=(Pindex.size + Qindex.size, tindex.size + vindex.size))
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
    print(actualpower)
    estpower=np.empty(0)
    for pi in range(0,Pindex.size):
        estpower=np.append(estpower,act_pow(Pindex[pi],v.size,v,teta,g,b))
    for qi in range(0,Qindex.size):
        estpower=np.append(estpower,react_pow(Qindex[qi],v.size,v,teta,g,b))
    print(estpower)
    return actualpower-estpower

def newtonrhapson(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b):
    jacobiinv=np.linalg.inv(jacobi(Pindex,Qindex,tindex,vindex,v,teta,g,b))
    deltapower= power_missmatch(Pactual,Qactual,Pindex,Qindex,vindex,v,teta,g,b)
    correction=jacobiinv.dot(deltapower)
    print('Jacobi-matrix:\n',jacobi(Pindex,Qindex,tindex,vindex,v,teta,g,b))
    print('\nVoltage angles:')
    for a in range(0,tindex.size):
        teta[a]=teta[a]+correction[a]
        print('Bus nr',tindex[a],': Voltage angle =',teta[a] )
    print('\nVoltage magnetudes:')
    for vm in range(0,vindex.size):
        v[vm-1]=v[vm]+correction[tindex.size-1+vm]
        print('Bus nr',vindex[vm],': Voltage magnitude =',v[vm])
    #for i in correction:
       # while i>0.5:
        #    newtonrhapson(Pactual,Qactual,Pindex,Qindex,tindex,vindex,v,teta,g,b)
    #return 0