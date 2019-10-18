import sys
sys.path.append(".")
from tests.basecase import *
from DCflow.distfactor import *
from tests.basecase import *

#print(dcflow(Pactual,Qactual,Pindex, Qindex, tindex, v, teta, g, b,z))
print(distfactor(Pindex,tindex,z,1,2))