import numpy as np
from Loadflow.gb import *
from Loadflow.tu import *

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
            row = np.append(row, -v[busi-1] * tij(busi, busj, teta, g, b))
    return row


def dqdt(busi, tindex, v, teta, g, b):
    row = np.empty(0)
    for busj in tindex:
        if busj==busi:
            element = 0
            for vj in range(1, v.size + 1):
                if vj != busi:
                    element += (
                            -v[busi - 1] * v[vj - 1] * tij(busi, vj, teta, g, b)
                    )
            row = np.append(row, element)

        else:
            row = np.append(
                row, v[busi - 1] * v[busj - 1] * tij(busi, busj, teta, g, b)
            )
    return row


def dqdv(busi, vindex, v, teta, g, b):
    row = np.empty(0)
    for busj in vindex:
        if busj == busi:
            element = 0
            for vj in range(1,v.size+1):
                if vj==busj:
                    element += -2 * v[busi-1] * b[busi-1][busi-1]
                else:
                    element += -v[vj-1] * uij(busi, vj, teta, g, b)
            row = np.append(row, element)
        else:
            row = np.append(row, -v[busi-1] * uij(busi, busj, teta, g, b))
    return row
