from tests.four_bus_system import *
from DCflow.dcflow import *
from Contingency.imml import *
#basecase
print('basecase: ')
pflow=dcflow(Pactual,Pindex, tindex, v, teta, g, b,z)
print(pflow)

#ex a)
print('line 23 is out:')
teta[:-1]=angles_imml(2,3,-4,Pactual,Pindex,tindex,z)
print('angles: ',teta)
print('flow 12:',line_flow(1,2,teta,z))
print('flow 13:',line_flow(1,3,teta,z))
print('flow 23:',line_flow(2,3,teta,z))
print('flow 43:',line_flow(4,3,teta,z))

#ex b)
#print('One of the lines 1-3 out:')
#teta[:-1]=angles_imml(1,3,10,Pactual,Pindex,tindex,z)
#print('angles: ',teta)
#print('flow 12:',line_flow(1,2,teta,z))
#print('flow 13:',line_flow(1,3,teta,z))
#print('flow 23:',line_flow(2,3,teta,z))
#print('flow 43:',line_flow(4,3,teta,z))