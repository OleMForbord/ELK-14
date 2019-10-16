import sys
sys.path.append(".")
from tests.basecase import *
from Loadflow.voltstab import *

#print(power_missmatch(Pactual,Qactual,Pindex,Qindex,v,teta,g,b))
#print_fdlf(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b,z,)
#print_primal(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b,z)
#print(h_matrix(Pindex,tindex,v,teta,g,b))
#print(b)
#newtonrhapson_print(Pactual,Qactual,Pindex,Qindex,tindex,vindex,v,teta,g,b)
#voltageStabilityFlatStart(Pactual,Qactual,Pindex, Qindex, tindex, vindex, g, b)
voltageStabilityFlatStart(Pactual,Qactual,Pindex, Qindex, tindex, vindex, g, b)


