import math
import sys
sys.path.append(".")
from Loadflow.voltstab import *


#CorrectorV_jac finds the jacobi with extension for the corrector-phase when the voltage is held constant
def correctorV_jac(Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta,busi):
    height = Pindex.size+ Qindex.size +1
    width = tindex.size + vindex.size +1
    halfHeight = math.floor(height/2)
    corrV = np.zeros(shape=(height,width))
    np.set_printoptions(suppress=True)
    corrV[0:height-1,0:width-1] = jacobi(Pindex, Qindex, tindex, vindex, v, teta, g, b)
    for i in range(height):
        if(i < halfHeight):
            corrV[i, width-1] = beta[i]
        elif(i<height-1):
            corrV[i, width-1] = alpha[i-halfHeight]
    corrV[height - 1, tindex.size + busi - 1] = 1
    return corrV

#CorrectorV_vector finds the correction of the angles when the voltage is held constant
def correctorV_vector(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta,busi):
    corrV_inv = np.linalg.inv(correctorV_jac(Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta,busi))
    missmatch = power_missmatch(Pactual,Qactual,Pindex,Qindex,v,teta,g,b)
    missmatch=np.append(missmatch,[0])
    corrV_vector=corrV_inv.dot(missmatch)
    return corrV_vector

#CorrectorV_phase uses the correctorV_jac and the correctorV_vector to change the values of the angels and vltages
def correctorV_phase(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta,step,busi):
    it = 0
    epsilonError = 0.001
    correction = correctorV_vector(Pactual, Qactual, Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta,busi)

    print('teta: ',teta)
    print('v: ',v)
    print('correction: ', correction)
    while np.any(abs(correction[:-1]) > epsilonError) == True:
        if (it >= 1000):  # CONVERGENCE_LIMIT):
            print("diverg")
            return 0
        for i in range(0, tindex.size):
            teta[tindex[i] - 1] += correction[i]
            print(teta)
        for j in range(0, vindex.size):
            v[vindex[j] - 1] += correction[tindex.size + j]
            print(v)
        correction = correctorV_vector(Pactual, Qactual, Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta,busi)
        it += 1
    return 1
