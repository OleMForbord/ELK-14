from Loadflow.PowerEq import *
import math



def predictor(Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta):
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

def correctorV(Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta):
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
        else:
            corrV[i,beta.size:width-1] = 1
    return corrV

def correctorP(Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta):
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

def predictor_vector(Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta):
    pred_inv = np.linalg.inv(predictor(Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta))
    missmatch=np.zeros(Pindex.size+Qindex.size)
    missmatch=np.append(missmatch,[1])
    predvector=pred_inv.dot(missmatch)
    predvector=np.delete(predvector,-1)
    return predvector

def correctorV_vector(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta):
    corrV_inv = np.linalg.inv(predictor(Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta))
    missmatch = power_missmatch(Pactual,Qactual,Pindex,Qindex,v,teta,g,b)
    missmatch=np.append(missmatch,[0])
    corrV_vector=corrV_inv.dot(missmatch)
    corrV_vector=np.delete(corrV_vector,-1)
    return corrV_vector

def correctorP_vector(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta):
    corrP_inv = np.linalg.inv(predictor(Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta))
    missmatch = power_missmatch(Pactual, Qactual, Pindex, Qindex, v, teta, g, b)
    missmatch = np.append(missmatch, [0])
    corrP_vector = corrP_inv.dot(missmatch)
    corrP_vector = np.delete(corrP_vector, -1)
    return corrP_vector

def contflow_print(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta):
    step=0.0786
    newtonrhapson(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b)
    print(predictor(Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha,beta))
    for busi in vindex:
        v[busi-1]+=





