import sys
sys.path.append(".")
import Loadflow.voltstab as lf

def missmat_p(Pactual,Qactual,Pindex,Qindex,tindex,v,teta,g,b):
    return lf.power_missmatch(Pactual,Qactual,Pindex,Qindex,v,teta,g,b)[:tindex.size]

def missmat_q(Pactual,Qactual,Pindex,Qindex,tindex,v,teta,g,b):
    return lf.power_missmatch(Pactual,Qactual,Pindex,Qindex,v,teta,g,b)[tindex.size:]
