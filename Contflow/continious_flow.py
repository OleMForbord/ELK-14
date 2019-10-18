import sys
sys.path.append(".")
import copy
from Contflow.predictor import *
from Contflow.corrector_s import *
from Contflow.corrector_v import *




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
            Pactual -= correction[-1] * beta
            Qactual -= correction[-1] * alpha
            print('\nCorrector phase, constant voltage\n')
            print(correctorV_jac(Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta,worstV_bus))
            print('Correction: ', correction)
            print('\nVoltage ang: ', teta)
            print('Voltage mag: ', v)
            sensitivity = predictor_vector(Pindex, Qindex,predictor_jac(Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta))
            if(check_sensitivity(sensitivity,tindex)==0):
               # teta += -teta + teta_init
               # v += -v +v_init
               # Pactual += -Pactual+Pactual_init
                #Qactual += -Qactual+Qactual_init

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








