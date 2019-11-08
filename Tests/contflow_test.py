import sys
sys.path.append(".")
from Contflow.continious_flow import *
from Cases.basecase import *

#newtonrhapson_print(Pactual,Qactual,Pindex,Qindex,tindex,vindex,v,teta,g,b)
contflow_print(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta,0.3)
#correctorV_phase(Pactual,Qactual,Pindex, Qindex, tindex, vindex, v, teta, g, b, alpha, beta,1,1)
#print(predictor_vector(Pindex,Qindex,predictor_jac(Pindex, Qindex, tindex, vindex, v, teta, g, b,alpha,beta)))
#print(check_sensitivity(predictor_vector(Pindex,Qindex,predictor_jac(Pindex, Qindex, tindex, vindex, v, teta, g, b,alpha,beta)),tindex))


