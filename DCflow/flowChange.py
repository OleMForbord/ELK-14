from DCflow.distfactor import *
#import DCflow.distfactor as dc
from numpy import *



def flowChange(Pactual, Pchange, Pindex, tindex, z, i,j):
    Pnew = Pactual - Pchange
    flowChange = 0
    for k in range (0,Pindex.size):
        flowChange += distfactor(Pindex, tindex, z, i,j)[k]*(Pnew[Pindex[k]-1])
    return flowChange
