import Loadflow.PowerEq as pe
import numpy as np


def h_matrix(Pindex,tindex,b):
    matrix= np.zeros(shape=(Pindex.size,tindex.size))
    for i in Pindex:
        row=np.zeros(tindex.size)
        for j in tindex:
            if(j!=i):
                row[j-1]=b[i-1][j-1]
            else:
                row[j - 1] = -b[i - 1][j - 1]
        matrix[i-1]=row
    return matrix

def l_matrix(Qindex,vindex,b):
    matrix= np.zeros(shape=(Qindex.size,vindex.size))
    for i in Qindex:
        row=np.zeros(vindex.size)
        for j in vindex:
            if(j!=i):
                row[j-1]=b[i-1][j-1]
            else:
                row[j - 1] = -b[i - 1][j - 1]
        matrix[i-1]=row
    return matrix

def heq_matrix(Pindex,tindex,z):
    matrix= np.zeros(shape=(Pindex.size,tindex.size))
    for i in Pindex:
        for j in tindex:
            if j!=i:
                matrix[i-1][j-1] = -1/(z[i-1][j-1].imag)
        for k in range(0, z[i-1].size):
            if z[i-1][k] != 0:
                matrix[i-1][i-1] += np.sum(1 / (z[i-1][k].imag))
    return matrix

def leq_matrix(Qindex,vindex,z):
    matrix= np.zeros(shape=(Qindex.size,vindex.size))
    for i in Qindex:
        for j in vindex:
            if j!=i:
                matrix[i-1][j-1] = -1/(z[i-1][j-1].imag)
        for k in range(0, z[i-1].size):
            if z[i-1][k] != 0:
                matrix[i-1][i-1] += np.sum(1 / (z[i-1][k].imag))
    return matrix

def missmat_p(Pactual,Qactual,Pindex,Qindex,tindex,v,teta,g,b):
    return pe.power_missmatch(Pactual,Qactual,Pindex,Qindex,v,teta,g,b)[:tindex.size]

def missmat_q(Pactual,Qactual,Pindex,Qindex,tindex,v,teta,g,b):
    return pe.power_missmatch(Pactual,Qactual,Pindex,Qindex,v,teta,g,b)[tindex.size:]

def print_standard(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b,z):
    it = 1
    epsilonError = 0.001
    matrix=np.zeros(shape=(Pindex.size+Qindex.size,tindex.size+vindex.size))
    matrix[:Pindex.size,0:tindex.size]=heq_matrix(Pindex,tindex,z)
    matrix[Pindex.size:,tindex.size:] = leq_matrix(Qindex,vindex,z)

    matrix_inv=np.linalg.inv(matrix)
    correction=matrix_inv.dot(pe.power_missmatch(Pactual,Qactual,Pindex,Qindex,v,teta,g,b))
    print('standard decoupled matrix:\n',matrix)
    print('\nVoltage angles: ', teta)
    print('Voltage magnetudes: ', v)

    while np.any(abs(correction)  > epsilonError):
        if(it>1000):
            print('the primal method did not converge, trying the dual:')
            return 0
        correction = matrix_inv.dot(pe.power_missmatch(Pactual, Qactual, Pindex, Qindex, v, teta, g,b))
        for a in range(0,tindex.size):
            teta[tindex[a]-1] += correction[a]
        for vm in range(0, vindex.size):
            v[vindex[vm]-1] += correction[tindex.size+vm]
        print('\niteraton: ', it)
        print('\nAngle correction:', correction[Pindex.size:])
        print('Voltage angles: ', teta)
        print('\nMagnitudes correction:',correction[:Pindex.size])
        print('Voltage magnetudes: ', v)




        correction= matrix_inv.dot(pe.power_missmatch(Pactual, Qactual, Pindex, Qindex, v, teta, g, b))

        it += 1
    it = 1
    return 1

def print_primal(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b,z):
    it = 1
    epsilonError = 0.001
    bm_inv=np.linalg.inv(bm(Pindex,tindex,b))
    bmm_inv=np.linalg.inv(bmm(Qindex,vindex,z))

    correction_teta = bm_inv.dot(missmat_p(Pactual, Qactual, Pindex, Qindex, tindex, v, teta, g, b))
    correction_v = bmm_inv.dot(missmat_q(Pactual, Qactual, Pindex, Qindex, vindex, v, teta, g, b))


    print('bm:\n',np.linalg.inv(bm_inv))
    print('\nbmm:\n',np.linalg.inv(bmm_inv))
    print('\nVoltage angles: ', teta)
    print('Voltage magnetudes: ', v)

    while np.any(abs(correction_teta)  > epsilonError):
        if(it>1000):
            print('the primal method did not converge, trying the dual:')
            return 0
        correction_teta = bm_inv.dot(missmat_p(Pactual, Qactual, Pindex, Qindex, tindex, v, teta, g, b))
        for a in range(0,tindex.size):
            teta[tindex[a]-1] += correction_teta[a]
        print('\niteraton: ', it)
        print('\nAngle correction:', correction_teta)
        print('Voltage angles: ', teta)

        correction_v = bmm_inv.dot(missmat_q(Pactual, Qactual, Pindex, Qindex, vindex, v, teta, g, b))
        for vm in range(0, vindex.size):
            v[vindex[vm]-1] += correction_v[vm]
        print('\nMagnitudes correction:',correction_v)
        print('Voltage magnetudes: ', v)

        it += 1
    it = 1
    return 1

def print_dual(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b,z):
    it = 1
    epsilonError = 0.001
    bm_inv=np.linalg.inv(bm(Pindex,tindex,b))
    bmm_inv=np.linalg.inv(bmm(Qindex,vindex,z))

    correction_teta = bm_inv.dot(missmat_p(Pactual, Qactual, Pindex, Qindex, tindex, v, teta, g, b))
    correction_v = bmm_inv.dot(missmat_q(Pactual, Qactual, Pindex, Qindex, vindex, v, teta, g, b))


    print('bm:\n',np.linalg.inv(bm_inv))
    print('\nbmm:\n',np.linalg.inv(bmm_inv))
    print('\nVoltage angles: ', teta)
    print('Voltage magnetudes: ', v)

    while np.any(abs(correction_teta)  > epsilonError):
        if it>1000:
            print('No solution found')
            return 0
        correction_v = bm_inv.dot(missmat_q(Pactual, Qactual, Pindex, Qindex, vindex, v, teta, g, b))
        for vm in range(0, vindex.size):
            v[vindex[vm] - 1] += correction_v[vm]
        print('\nMagnitudes correction:', correction_v)
        print('Voltage magnetudes: ', v)

        correction_teta = bmm_inv.dot(missmat_p(Pactual, Qactual, Pindex, Qindex, tindex, v, teta, g, b))
        for a in range(0,tindex.size):
            teta[tindex[a]-1] += correction_teta[a]
        print('\niteraton: ', it)
        print('\nAngle correction:', correction_teta)
        print('Voltage angles: ', teta)
        it+=1
    return 1

def print_fdlf(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b,z):
    if print_primal(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b,z)==0:
        print_dual(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b,z)

def dcflow(Pactual,Qactual,Pindex, Qindex, tindex, v, teta, g, b,z):
    print('initial angles:\n',teta)
    b_matrix=heq_matrix(Pindex,tindex,z)
    print('B-matrix:\n',b_matrix)
    bmm_inv=np.linalg.inv(b_matrix)
    correction=bmm_inv.dot(missmat_p(Pactual,Qactual,Pindex,Qindex,tindex,v,teta,g,b))
    print('correction: ',correction)
    for a in range(0, tindex.size):
        teta[tindex[a] - 1] += correction[a]
    print('Updated angles: ',teta)
    pflow=np.zeros(shape=(teta.size,teta.size))
    for i in range(0,teta.size):
        for j in range(0,teta.size):
            if z[i][j]!=0:
                pflow[i][j]=(teta[i]-teta[j])/z[i][j].imag
    return pflow

def distfactor(Pindex,tindex,b,z,i,j):

    heq=heq_matrix(Pindex,tindex,z)
    heq_new = np.array([heq[i-1],heq[j-1]])
    print(heq_new)
    return np.linalg.inv(heq_new).dot([1/z[i-1][j-1].imag,-1/z[i-1][j-1].imag])