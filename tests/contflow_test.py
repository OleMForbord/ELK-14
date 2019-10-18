import sys
sys.path.append(".")
from Contflow.continious_flow import *
from tests.basecase import *

#newtonrhapson_print(Pactual,Qactual,Pindex,Qindex,tindex,vindex,v,teta,g,b)
contflow_print(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta,0.2)
#correctorV_phase(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta,1,1)

