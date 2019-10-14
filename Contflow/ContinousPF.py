from Loadflow.PowerEq import *
import math
import copy

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
    #corrV_vector=np.delete(corrV_vector,-1)
    return corrV_vector

#CorrectorV_phase uses the correctorV_jac and the correctorV_vector to change the values of the angels and vltages
def correctorV_phase(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta,step,busi):
    it = 0
    epsilonError = 0.001
    correction = correctorV_vector(Pactual, Qactual, Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta,busi)
    while np.any(abs(correction[:-1]) > epsilonError) == True:
        if (it >= 1000):  # CONVERGENCE_LIMIT):
            print("diverg")
            return 0
        for i in range(0, tindex.size):
            teta[tindex[i] - 1] += correction[i]
        for j in range(0, vindex.size):
            v[vindex[j] - 1] += correction[tindex.size + j]
        Pactual += step * beta
        Qactual += step * alpha
        correction = correctorV_vector(Pactual, Qactual, Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta,busi)
        it += 1
    return 1

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
    teta_init=copy.copy(teta)
    v_init=copy.copy(v)
    #teta_init=np.zeros(teta.size)
    #teta_init+=teta
    #v_init=np.zeros(v.size)
    #v_init+=v
    it = 0
    epsilonError = 0.001
    correction = correctorS_vector(Pactual, Qactual, Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta)
    while np.any(abs(correction) > epsilonError) == True:
        if (it >= 1000):  # CONVERGENCE_LIMIT):
            teta+=-teta+teta_init
            v+=-v+v_init
            print("diverg")
            return 0
        for i in range(0, tindex.size):
            teta[tindex[i] - 1] += correction[i]
        for j in range(0, vindex.size):
            v[vindex[j] - 1] += correction[tindex.size + j]
        correction = correctorS_vector(Pactual, Qactual, Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta)
        it += 1
    return 1

#Check sensitivity returns TRUE if all the voltage-sensitivities are negative or zero, 0 if not.
def check_sensitivity(sensitivity,tindex):
    return np.all(sensitivity[tindex.size:]<=0)


def contflow_print(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta,step):
    it=1
    newtonrhapson(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b)
    P_load = np.array(sum(abs(Pactual)))
    v1 = np.array([v[0]])
    v2 = np.array([v[1]])
    sensitivity=predictor_vector(Pindex,Qindex,predictor_jac(Pindex, Qindex, tindex, vindex, v, teta, g, b,alpha,beta))
    while (check_sensitivity(sensitivity,tindex)):
        Pactual_init=copy.copy(Pactual)
        Qactual_init=copy.copy(Qactual)
        teta_init = copy.copy(teta)
        v_init=copy.copy(v)


        print('\niteration: ',it)
        print('\nVoltage ang: ', teta)
        print('Voltage mag: ', v)
        print('\nActive power injections: ', pinj(v, teta, g, b))
        print('Reactive power injection: ', qinj(v, teta, g, b))

        pred_jac = predictor_jac(Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta)
        print('\nPredictor phase:\n')
        print(pred_jac)
        print('sensitivity: ', predictor_vector(Pindex, Qindex, pred_jac))

        pred_phase(Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta, step)
        print('\nVoltage angles est: ', teta)
        print('Voltage mag    est: ', v)

        Pactual -= step * beta
        Qactual -= step * alpha

        if(correctorS_phase(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b,alpha,beta)):
            correction=correctorS_vector(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b,alpha,beta)
            print('\nCorrector phase, constant load\n')
            print(correctorS_jac(Pindex, Qindex, tindex, vindex, v, teta, g, b,alpha,beta))
            print('Correction: ',correction)
            print('\nVoltage ang: ', teta)
            print('Voltage mag: ', v)
        else:
            worstV_bus=np.argmax(sensitivity[tindex.size:])+1

            correction = correctorV_vector(Pactual, Qactual, Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta,worstV_bus)
            correctorV_phase(Pactual, Qactual, Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta,step,worstV_bus)

            print('\nCorrector phase, constant voltage\n')
            print(correctorV_jac(Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta,worstV_bus))
            print('Correction: ', correction)
            print('\nVoltage ang: ', teta)
            print('Voltage mag: ', v)
            sensitivity = predictor_vector(Pindex, Qindex,predictor_jac(Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta))
            if(check_sensitivity(sensitivity,tindex)==0):
                teta += -teta + teta_init
                v += -v +v_init
                Pactual += -Pactual+Pactual_init
                Qactual += -Qactual+Qactual_init

                print('\nPrevious valid solution:')
                print('Voltage ang: ', teta)
                print('Voltage mag: ', v)
                P_load = np.append(P_load, sum(abs(Pactual)))
                v1 = np.append(v1, v[0])
                v2 = np.append(v2, v[1])
                break

        sensitivity = predictor_vector(Pindex, Qindex, predictor_jac(Pindex, Qindex, tindex, vindex, v, teta, g, b,alpha,beta))
        P_load = np.append(P_load,sum(abs(Pactual)))
        v1 = np.append(v1, v[0])
        v2 = np.append(v2, v[1])
        print('sensitivity: ',sensitivity)
        it+=1
    plt.xlabel("Load power")
    plt.ylabel("Voltage")
    plt.plot(P_load, v1)
    plt.plot(P_load, v2)
    plt.show()








