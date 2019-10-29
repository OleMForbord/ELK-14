from tests.four_bus_system import *
from DCflow.dcflow import *
from Contingency.imml import *
#basecase
print('basecase: ')
pflow=dcflow(Pactual,Pindex, tindex, v, teta, g, b,z)
print(pflow)

#ex a)
#print('line 23 is out:')
#lineout = [[2,3]] #lineout is a list of lines completely out
#teta[:-1]=angles_imml(2,3,-4,Pactual,Pindex,tindex,z)
#print('angles: ',teta)
#print_systemflow(teta,z,lineout)
#ex b)
lineout=[]
print_systemflow(teta, z, lineout)

