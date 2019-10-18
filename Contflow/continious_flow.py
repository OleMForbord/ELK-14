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
    print(sensitivity)
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
            worstV_bus=vindex[np.argmax(abs(sensitivity[tindex.size:]))]
            print('worst bus: ',worstV_bus)

            correction = correctorV_vector(Pactual, Qactual, Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta,worstV_bus)
            correctorV_phase(Pactual, Qactual, Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta,step,worstV_bus)
            print('Pactual before corr:',Pactual)
            print('Qactual before corr:', Qactual)
            Pactual -= correction[-1] * beta
            Qactual -= correction[-1] * alpha
            print('Pactual before corr:', Pactual)
            print('Qactual before corr:', Qactual)
            print('\nCorrector phase, constant voltage\n')
            print(correctorV_jac(Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta,worstV_bus))
            print('Correction: ', correction)
            print('\nVoltage ang: ', teta)
            print('Voltage mag: ', v)
            
            sensitivity = predictor_vector(Pindex, Qindex,predictor_jac(Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta))
            if(check_sensitivity(sensitivity,tindex)==0 or np.any(v<=0)):
                teta += teta_init-teta
                v += v_init-v-teta
                Pactual += Pactual_init-Pactual
                Qactual += Qactual_init-Qactual
                v1=v1[:-1]
                v2=v2[:-1]
                P_load=P_load[:-1]
                print('\nPrevious valid solution:')
                print('Voltage ang: ', teta)
                print('Voltage mag: ', v)
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
    plt.show()







