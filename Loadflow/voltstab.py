import sys
sys.path.append(".")

from Loadflow.newtonrapshon import *

def voltageStabilityFlatStart(Pactual,Qactual,Pindex, Qindex, tindex, vindex, g, b):
    syst_load=abs(Pactual[0]+Pactual[1])
    P_load=np.zeros(0)
    v1 = np.zeros(0)
    v2 = np.zeros(0)

    while True:
        teta = np.array([0.0, 0.0, 0.0])  # creates an array with the initial angles
        v = np.array([1.0, 1.0, 1.0])  # creates an array with the initial voltages

        if (newtonrhapson(Pactual, Qactual, Pindex, Qindex, tindex, vindex, v, teta, g, b) == 0):
            print("breakes")
            break
        Pactual = Pactual - [0.2 * 0.3, 0.2 * 0.7]
        # Qactual = Qactual - [0.2 * 0.3, 0.2 * 0.7]
        syst_load += 0.2
        print(syst_load)
        P_load = np.append(P_load, syst_load)
        v1 = np.append(v1, v[0])
        v2 = np.append(v2, v[1])
    print("Divergence at Pactual: ", Pactual)
    print(v1, v2, P_load)

    plt.xlabel("Load power")
    plt.ylabel("Voltage")
    plt.plot(P_load, v1)
    plt.plot(P_load, v2)
    plt.show()

def voltageStabilityAccumulating(Pactual,Qactual,Pindex, Qindex, tindex, vindex, teta, v,  g, b):
    while True:
        print("Teta is:", teta)
        if (newtonrhapson(Pactual, Qactual, Pindex, Qindex, tindex, vindex, v, teta, g, b) == 0):
            break
            Pactual = Pactual - [0.2 * 0.3, 0.2 * 0.7]
        # Qactual = Qactual - [0.2 * 0.3, 0.2 * 0.7]

    print("Divergence at Pactual: ", Pactual)
    plt.xlabel("Load power")
    plt.ylabel("Voltage")
    plt.show()