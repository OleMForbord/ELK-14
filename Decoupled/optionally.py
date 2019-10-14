import Loadflow.PowerEq as pe
import numpy as np


def h_matrix(Pindex, tindex, v, teta, g, b):
    matrix=np.zeros(shape=(Pindex.size,tindex.size))
    for i in Pindex:
            matrix[i-1]=pe.dpdt(i,tindex,v,teta,g,b)
    return matrix

def n_matrix(Pindex, vindex, v, teta, g, b):
    matrix = np.zeros(shape=(Pindex.size, vindex.size))
    for i in Pindex:
        matrix[i - 1] = pe.dpdv(i, vindex, v, teta, g, b)
    return matrix

def m_matrix(Qindex, tindex, v, teta, g, b):
    matrix = np.zeros(shape=(Qindex.size, tindex.size))
    for i in Qindex:
        matrix[i - 1] = pe.dqdt(i, tindex, v, teta, g, b)
    return matrix

def l_matrix(Qindex, vindex, v, teta, g, b):
    matrix = np.zeros(shape=(Qindex.size, vindex.size))
    for i in Qindex:
        matrix[i - 1] = pe.dqdv(i, vindex, v, teta, g, b)
    return matrix

def heq_matrix(Pindex,Qindex,tindex,vindex,v,teta,g,b):
    n=n_matrix(Pindex, vindex, v, teta, g, b)
    l_inv=np.linalg.inv(l_matrix(Qindex, vindex, v, teta, g, b))
    m=m_matrix(Qindex, tindex, v, teta, g, b)
    return h_matrix(Pindex,tindex,v,teta,g,b)-n.dot(l_inv).dot(m)

def leq_matrix(Pindex,Qindex,tindex,vindex,v,teta,g,b):
    m = m_matrix(Qindex, tindex, v, teta, g, b)
    h_inv=np.linalg.inv(h_matrix(Pindex, tindex, v, teta, g, b))
    n = n_matrix(Pindex, vindex, v, teta, g, b)
    return l_matrix(Pindex, vindex, v, teta, g, b) - m.dot(h_inv).dot(n)