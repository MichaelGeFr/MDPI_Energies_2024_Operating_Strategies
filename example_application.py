#!/usr/bin/env python
# -*-coding:utf-8 -*-

'''
__author__      :   Lukas Theisinger 
__version__     :   1.0
__maintainer__  :   Michael Frank
__contact__     :   m.frank@ptw.tu-darmstadt.de
'''


from pyomo.environ import *
from RuBOS.optimization.ControlLogic import *

cl_data = {
    None: {
        'I': {None: [1]}, # influence factor index
        'J': {None: [1, 2]}, # decision rule index
        'T': {None: [0, 1]}, # time index
        'P': {None: [1, 2]}, # prio index
        'C': {None: [44, 88]}, # converter index
        'A': {1: 1},
        'B': {1: 1, 2: 2},
        'P_nom': {44: 10, 88: 10},
    }
}

model_data = {
    None: {
        'I': {None: [1]}, # influence factor index
        'J': {None: [1, 2]}, # decision rule index
        'T': {None: [0, 1]}, # time index
        'P': {None: [1, 2]}, # prio index
        'C': {None: [44, 88]}, # converter index
        'A': {1: 1},
        'B': {1: 1, 2: 2},
        'P_nom': {44: 10, 88: 10},
        'P_dem': {0: 12, 1: 5},
        'c_el': {None: 30},
        'eta': {(44, 0): 5, (44, 1): 2, (88, 0): 3, (88, 1): 3}
    }
}

model = AbstractModel()

### indices ###
model.I = Set(domain = NonNegativeIntegers, doc = "influence factor index I")
model.T = Set(domain = NonNegativeIntegers, doc = "time step index T")
model.C = Set(doc = "converter index C")

### variables ###
model.P_gen = Var(model.C, model.T)

### parameters ###
model.P_nom = Param(model.C)
model.P_dem = Param(model.T)
model.eta = Param(model.C, model.T)
model.c_el = Param()

model.controlLogic = ControlLogic.create_instance(cl_data)

def influence_factor_def(m, I, T):
    return m.controlLogic.fac[I, T] == m.eta[44, T]
model.infl_fac_def = Constraint(model.I, model.T, rule = influence_factor_def)

def p_gen_connection(m, C, T):
    return m.controlLogic.P_gen[C, T] == m.P_gen[C, T]
model.p_gen_con = Constraint(model.C, model.T, rule = p_gen_connection)

def power_limit(m, C, T):
    return m.P_gen[C, T] <= m.P_nom[C]
model.power_limit_cons = Constraint(model.C, model.T, rule = power_limit)

def energy_balance(m, T):
    return m.P_dem[T] == sum(m.P_gen[c, T] for c in m.C)
model.energy_balance_cons = Constraint(model.T, rule = energy_balance)

def obj(m):
    return sum((m.P_gen[c, t] * m.c_el) / m.eta[c, t] for t in m.T for c in m.C)
model.obj_rule = Objective(rule = obj)

opt = SolverFactory('cplex')
instance = model.create_instance(model_data)
opt.solve(instance, tee = True)
instance.display()








