from PowerEq import *
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
