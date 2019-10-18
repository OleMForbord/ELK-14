import math
import sys
sys.path.append(".")
from Loadflow.voltstab import *

#CorrectorS_jac finds the jacobi with extension for the corrector-phase when the load is held constant
def correctorS_jac(Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta):
    height = Pindex.size+ Qindex.size +1
    width = tindex.size + vindex.size +1
    halfHeight = math.floor(height/2)
    corrP = np.zeros(shape=(height,width))
    np.set_printoptions(suppress=True)
    corrP[0:height-1,0:width-1] = jacobi(Pindex, Qindex, tindex, vindex, v, teta, g, b)

    for i in range(height):
        if(i < halfHeight):
            corrP[i, width-1] = beta[i]
        elif(i<height-1):
            corrP[i, width-1] = alpha[i-halfHeight]
        else:
            corrP[i,width-1] = 1
    return corrP

#CorrectorS_vector finds the correctio of the angles and voltages, when load is held constant
def correctorS_vector(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta):
    corrP_inv = np.linalg.inv(correctorS_jac(Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta))
    missmatch = power_missmatch(Pactual, Qactual, Pindex, Qindex, v, teta, g, b)
    missmatch = np.append(missmatch, [0])
    corrP_vector = corrP_inv.dot(missmatch)
    corrP_vector = np.delete(corrP_vector, -1)
    return corrP_vector

#CorrectorS_phase uses the correctorS_jac and the correctorS_vector to change the values of the angels and vpltages
def correctorS_phase(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta):
    #teta_init=copy.copy(teta)
    #v_init=copy.copy(v)
    teta_init=np.zeros(teta.size)
    teta_init+=teta
    v_init=np.zeros(v.size)
    v_init+=v
    it = 0
    epsilonError = 0.001
    correction = correctorS_vector(Pactual, Qactual, Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta)
    while np.any(abs(correction) > epsilonError) == True:
        if (it >= 1000):  # CONVERGENCE_LIMIT):
            teta+=teta_init-teta
            v+=v_init-v
            print("diverg")
            print('teta:',teta)
            print('v',v)
            return 0
        for i in range(0, tindex.size):
            teta[tindex[i] - 1] += correction[i]
        for j in range(0, vindex.size):
            v[vindex[j] - 1] += correction[tindex.size + j]
        correction = correctorS_vector(Pactual, Qactual, Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta)
        it += 1
    return 1
