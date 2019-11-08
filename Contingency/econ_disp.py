import sys
sys.path.append(".")
from tests.four_bus_system import *
from DCflow.distfactor import *
import time
import pandas as pd
import pyomo.environ as pyo
from pyomo.opt import SolverFactory

'''
Some pyomo documentation:
Objective: expression to be max/minimized
Variables: represent the unkown part of the model, what we want to optimize
Parameters: data that must be suppliied to solve the problem
Constraint: constraint expressions of the model

steps of using pyomo as optimizer:
1) Create model and declare components
2) Instantiate the model
3) Apply solver
4) Interrogate solver results

To use glpk you had to add the pat of the executables(winglpk-4.65\glpk-4.65\w64) to environment variables
To use gerobi add path to the gerobi811/win64/bin
'''

'''
Gen=["G1","G2","G3","G4"] #Creates a list of differetn generators
Cost={"G1":1,"G2":10,"G3":2,"G4":3} #Creates a dictionary with key:values

ptdf_34={"G1":distfactor(Pindex,tindex,z,3,4),"G2":distfactor(Pindex,tindex,z,3,4),"G3":distfactor(Pindex,tindex,z,3,4),"G4":distfactor(Pindex,tindex,z,3,4)}

'''
model= pyo.ConcreteModel() #Create a model variable

model.G1=pyo.Var(within = pyo.NonNegativeReals)
model.G2=pyo.Var(within = pyo.NonNegativeReals)
model.G3=pyo.Var(within = pyo.NonNegativeReals)
model.G4=pyo.Var(within = pyo.NonNegativeReals)
ptdf34=np.array(distfactor(Pindex,tindex,z,3,4,slackbusnr))

c=np.array([3,10,2,1])

def Objective(model):
    #returns the objective function
    return(3*model.G1+10*model.G2+2*model.G3+model.G4)


model.OBJ = pyo.Objective(rule= Objective, sense= pyo.minimize)

def power_balance(model):
    return (model.G1+model.G2+model.G3+model.G4==np.sum(Pactual))

model.C1=pyo.Constraint(rule=power_balance)

def maxflow_34(model):
    return (model.G1*ptdf34[0]+model.G2*ptdf34[1]+model.G3*ptdf34[2]+model.G4*ptdf34[3]<=1.0)
model.C2=pyo.Constraint(rule=maxflow_34)

def minflow_34(model):
    return (-1.0<=model.G1*ptdf34[0]+model.G2*ptdf34[1]+model.G3*ptdf34[2]+model.G4*ptdf34[3])
model.C3=pyo.Constraint(rule=minflow_34)

opt = pyo.SolverFactory("glpk")
model.dual=pyo.Suffix(direction=pyo.Suffix.IMPORT)
result=opt.solve(model,load_solutions=True)

model.display()
model.dual.display()

