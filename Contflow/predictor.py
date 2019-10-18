import math
import sys
sys.path.append(".")
from Loadflow.voltstab import *

#Predictor_jac finds the jacobi with extension for the predicor-phase
def predictor_jac(Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta):
    height = Pindex.size+ Qindex.size +1
    width = tindex.size + vindex.size +1
    halfHeight = math.floor(height/2)
    pred = np.zeros(shape=(height,width))
    np.set_printoptions(suppress=True)
    pred[0:height-1,0:width-1] = jacobi(Pindex, Qindex, tindex, vindex, v, teta, g, b)
    for i in range(height):
        if(i < halfHeight):
            pred[i, width-1] = beta[i]
        elif(i<height-1):
            pred[i, width-1] = alpha[i-halfHeight]
        else:
            pred[i,width-1] = 1
    return pred

#Pred_vector finds the gradient/sensitivity of 1pu change in the load
def predictor_vector(Pindex, Qindex, pred_jac):
    pred_inv = np.linalg.inv(pred_jac)
    missmatch=np.zeros(Pindex.size+Qindex.size)
    missmatch=np.append(missmatch,[1])
    predvector=np.dot(pred_inv,missmatch)
    predvector=np.delete(predvector,-1)
    return predvector

#predictor_phase uses the gradient to change the values of teta and v, a step must be chosen
def pred_phase(Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta,step):
    pred_jac = predictor_jac(Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta)
    sensitivity = predictor_vector(Pindex, Qindex, pred_jac)
    for i in range(0,tindex.size):
        teta[tindex[i]-1] += step * sensitivity[i]
    for j in range(0,vindex.size):
        v[vindex[j]-1] += step * sensitivity[tindex.size+j]
    return sensitivity
