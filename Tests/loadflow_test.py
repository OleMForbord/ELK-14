import sys
sys.path.append(".")
from Cases.basecase import *
from Loadflow.voltstab import *
import cProfile

print(g_matrix(z))
print(b_matrix(z))
#print(power_missmatch(Pactual,Qactual,Pindex,Qindex,v,teta,g,b))
newtonrhapson_print(Pactual,Qactual,Pindex,Qindex,tindex,vindex,v,teta,g,b)
#voltageStabilityFlatStart(Pactual,Qactual,Pindex, Qindex, tindex, vindex, g, b)
#voltageStabilityFlatStart(Pactual,Qactual,Pindex, Qindex, tindex, vindex, g, b)


cProfile.run('''newtonrhapson_print(Pactual,Qactual,Pindex,Qindex,tindex,vindex,v,teta,g,b)''')

