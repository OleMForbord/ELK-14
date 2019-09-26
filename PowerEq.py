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
    estpower=np.empty(0)
    for pi in range(0,Pindex.size):
        estpower=np.append(estpower,act_pow(Pindex[pi],v.size,v,teta,g,b))
    print('active power injection: ', estpower)
    for qi in range(0,Qindex.size):
        estpower=np.append(estpower,react_pow(Qindex[qi],v.size,v,teta,g,b))
    print('reactive power injection: ',estpower[:Pindex.size])
    return actualpower-estpower

def newtonrhapson(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b,error_size):
    print('\nVoltage angles:', teta)
    print('Voltage magnitudes: ', v)
    jacobiinv=np.linalg.inv(jacobi(Pindex,Qindex,tindex,vindex,v,teta,g,b))
    missmatch= power_missmatch(Pactual,Qactual,Pindex,Qindex,v,teta,g,b)
    correction=jacobiinv.dot(missmatch)
    print('\nActive power missmatch: ',missmatch[:Pindex.size])
    print('Reactive power missmatch: ', missmatch[Pindex.size:])
    print('\nJacobi-matrix:\n',jacobi(Pindex,Qindex,tindex,vindex,v,teta,g,b))
    print('\ncorrection: ',correction)
    for a in range(0,tindex.size):
        teta[a]=teta[a]+correction[a]
    for vm in range(0,vindex.size):
        v[vm]=v[vm]+correction[tindex.size+vm]
    for i in range(0,3):
        if(abs(correction[1])>error_size):
            newtonrhapson(Pactual,Qactual,Pindex,Qindex,tindex,vindex,v,teta,g,b,error_size)
