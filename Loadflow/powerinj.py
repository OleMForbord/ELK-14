import sys
sys.path.append(".")
from Loadflow.direvatives import *

def act_pow(busi, n, v, teta, g, b):  # n=number of buses
    p = (v[busi - 1] ** 2) * g[busi - 1][busi - 1]
    for busj in range(1, n + 1):
        if busj == busi:
            continue
        else:
            p = p - v[busi - 1] * v[busj - 1] * tij(busi, busj, teta, g, b)
    return p


def react_pow(busi, n, v, teta, g, b):  # n=number of buses
    q = -(v[busi - 1] ** 2) * b[busi - 1][busi - 1]
    for busj in range(1, n + 1):
        if busj == busi:
            continue
        else:
            q = q - v[busi - 1] * v[busj - 1] * uij(busi, busj, teta, g, b)
    return q

def pinj(v,teta,g,b):
    powerinj=np.zeros(0)
    for busj in range(1,v.size+1):
        powerinj=np.append(powerinj,act_pow(busj,v.size,v,teta,g,b))
    return powerinj

def qinj(v,teta,g,b):
    powerinj=np.zeros(0)
    for busj in range(1,v.size+1):
        powerinj=np.append(powerinj,react_pow(busj,v.size,v,teta,g,b))
    return powerinj
def netinj(v,teta,g,b):
    return np.append(pinj(v,teta,g,b),qinj(v,teta,g,b))