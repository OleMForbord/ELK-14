import sys
sys.path.append(".")
from DCflow.distfactor import *
from DCflow.dcflow import *
from tests.basecase import *

print(dcflow(Pactual,Pindex,tindex,v,teta,g,b,z))
#print(distfactor(Pindex,tindex,z,1,2))