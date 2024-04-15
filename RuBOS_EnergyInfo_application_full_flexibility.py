#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
__author__      :   Lukas Theisinger 
__version__     :   1.0
__maintainer__  :   Michael Frank
__contact__     :   m.frank@ptw.tu-darmstadt.de
'''
import os
from pyomo.environ import *
from RuBOS.optimization.ControlLogic_validation import *
from RuBOS.utilities.utilities import *

current_dir=os.getcwd()

scen = "summer"
type = "flex"

indexing_dict = {
    'T': None,
    'c_el': 'T',
    'P_th_heat': 'T',
    'P_th_cool': 'T',
    'T_amb': 'T',
    'T_amb_avg': 'T',
    'T_flow_heat': 'T',
    'T_flow_cool': 'T'
}

model_data = prepare_timeseries(os.path.join(current_dir, "MDPI_data_Industry_" + scen + ".xlsx"), indexing_dict=indexing_dict)


""" 
data: P_th_heat, P_th_cool, T_amb, T_amb_avg, T_flow_heat, T_flow_cool, c_el

"""

model = AbstractModel()

### indices ###
model.T = Set(domain = NonNegativeIntegers, doc = "time step index T")
model.C = Set(doc = "converter index C")

### variables ###
model.P_el_hp1 = Var(model.T, bounds = (0, 375))
model.P_el_hp2 = Var(model.T, bounds = (0, 250))
model.P_el_hp3 = Var(model.T, bounds = (0, 125))
model.P_th_heat_hp1 = Var(model.T, within = NonNegativeReals)
model.P_th_heat_hp2 = Var(model.T, within = NonNegativeReals)
model.P_th_heat_hp3 = Var(model.T, within = NonNegativeReals)
model.P_th_cool_hp1 = Var(model.T, within = NonNegativeReals)
model.P_th_cool_hp2 = Var(model.T, within = NonNegativeReals)
model.P_th_cool_hp3 = Var(model.T, within = NonNegativeReals)
model.bOn_hp1 = Var(model.T, within = Binary)
model.bOn_hp2 = Var(model.T, within = Binary)
model.bOn_hp3 = Var(model.T, within = Binary)
model.bSwitchOn_hp1 = Var(model.T, within = Binary)
model.bSwitchOn_hp2 = Var(model.T, within = Binary)
model.bSwitchOn_hp3 = Var(model.T, within = Binary)

model.P_gas_chp1 = Var(model.T, bounds = (0, 6000))
model.P_el_chp1 = Var(model.T, within = NonNegativeReals)
model.P_th_heat_chp1 = Var(model.T, within = NonNegativeReals)
model.bOn_chp1 = Var(model.T, within = Binary)
model.bSwitchOn_chp1 = Var(model.T, within = Binary)

# model.P_gas_chp2 = Var(model.T, bounds = (0, 3000))
# model.P_el_chp2 = Var(model.T, within = NonNegativeReals)
# model.P_th_heat_chp2 = Var(model.T, within = NonNegativeReals)
# model.bOn_chp2 = Var(model.T, within = Binary)
# model.bSwitchOn_chp2 = Var(model.T, within = Binary)

model.P_el_ct = Var(model.T, bounds = (0, 10000))
model.P_th_cool_ct = Var(model.T, within = NonNegativeReals)

model.E_th_heat = Var(model.T, bounds = (0, 800))
model.E_th_cool = Var(model.T, bounds = (0, 800))

### parameters ###
model.corr_fac_hp = 0.5
model.rel_min_hp = 0.5
model.switch_loss_hp = 0.2

model.eta_th_chp = 0.5
model.eta_el_chp = 0.4
model.rel_min_chp = 0.5
model.switch_loss_chp = 0.2

model.T_amb = Param(model.T)
model.T_amb_avg = Param(model.T)
model.T_flow_heat = Param(model.T) # Kelvin
model.T_flow_cool = Param(model.T) # Kelvin
model.c_el = Param(model.T)
model.c_gas = 0.08 # €/kWh

model.P_th_heat = Param(model.T)
model.P_th_cool = Param(model.T)
model.E_th_heat_target = 800
model.E_th_cool_target = 800
model.dT = 1

def heat_balance1(m, T):
    if T == m.T.first():
        return m.E_th_heat[T] == m.E_th_heat_target + m.P_th_heat_hp1[T] + m.P_th_heat_hp2[T] + m.P_th_heat_hp3[T] + m.P_th_heat_chp1[T] - m.P_th_heat[T]
    else:
        return Constraint.Skip
model.heat_bal1_cons = Constraint(model.T, rule = heat_balance1)

def heat_balance2(m, T):
    if not T == m.T.first():
        return m.E_th_heat[T] - m.E_th_heat[T - 1] == m.P_th_heat_hp1[T] + m.P_th_heat_hp2[T] + m.P_th_heat_hp3[T] + m.P_th_heat_chp1[T] - m.P_th_heat[T]
    else:
        return Constraint.Skip
model.heat_bal2_cons = Constraint(model.T, rule = heat_balance2)

def cool_balance1(m, T):
    if T == m.T.first():
        return m.E_th_cool[T] == m.E_th_cool_target + m.P_th_cool_hp1[T] + m.P_th_cool_hp2[T] + m.P_th_cool_hp3[T] + m.P_th_cool_ct[T] - m.P_th_cool[T]
    else:
        return Constraint.Skip
model.cool_bal1_cons = Constraint(model.T, rule = cool_balance1)

def cool_balance2(m, T):
    if not T == m.T.first():
        return m.E_th_cool[T] - m.E_th_cool[T] == m.P_th_cool_hp1[T] + m.P_th_cool_hp2[T] + m.P_th_cool_hp3[T] + m.P_th_cool_ct[T] - m.P_th_cool[T]
    else:
        return Constraint.Skip
model.cool_bal2_cons = Constraint(model.T, rule = cool_balance2)

def hp1_balance1(m, T):
    return m.P_th_cool_hp1[T] + m.P_el_hp1[T] == m.P_th_heat_hp1[T]
model.hp1_bal_cons1 = Constraint(model.T, rule = hp1_balance1)

def hp1_balance2(m, T):
    return ((m.T_flow_heat[T] * m.corr_fac_hp)/(m.T_flow_heat[T] - m.T_flow_cool[T])) * m.P_el_hp1[T] - m.bSwitchOn_hp1[T] * m.switch_loss_hp * 375 == m.P_th_heat_hp1[T]
model.hp1_balance2 = Constraint(model.T, rule = hp1_balance2)

def hp1_balance3(m, T):
    return m.bOn_hp1[T] * 375 >= m.P_el_hp1[T]
model.hp1_bal_cons3 = Constraint(model.T, rule = hp1_balance3)

def hp1_balance4(m, T):
    return m.P_el_hp1[T] >= m.bOn_hp1[T] * m.rel_min_hp * 375
model.hp1_bal_cons4 = Constraint(model.T, rule = hp1_balance4)

def hp1_balance5(m, T):
    if not T == m.T.first():
        return m.bSwitchOn_hp1[T] >= m.bOn_hp1[T] - m.bOn_hp1[T-1]
    else:
        return Constraint.Skip
model.hp1_bal_cons5 = Constraint(model.T, rule = hp1_balance5)

def hp1_balance6(m, T):
    if not T == m.T.first():
        return m.bSwitchOn_hp1[T] + m.bSwitchOn_hp1[T - 1] <= 1
    else:
        return Constraint.Skip
model.hp1_bal_cons6 = Constraint(model.T, rule = hp1_balance6)

def hp2_balance1(m, T):
    return m.P_th_cool_hp2[T] + m.P_el_hp2[T] == m.P_th_heat_hp2[T]
model.hp2_bal_cons1 = Constraint(model.T, rule = hp2_balance1)

def hp2_balance2(m, T):
    return ((m.T_flow_heat[T] * m.corr_fac_hp)/(m.T_flow_heat[T] - m.T_flow_cool[T])) * m.P_el_hp2[T] - m.bSwitchOn_hp2[T] * m.switch_loss_hp * 250 == m.P_th_heat_hp2[T] 
model.hp2_balance2 = Constraint(model.T, rule = hp2_balance2)

def hp2_balance3(m, T):
    return m.bOn_hp2[T] * 250 >= m.P_el_hp2[T]
model.hp2_bal_cons3 = Constraint(model.T, rule = hp2_balance3)

def hp2_balance4(m, T):
    return m.P_el_hp2[T] >= m.bOn_hp2[T] * m.rel_min_hp * 250
model.hp2_bal_cons4 = Constraint(model.T, rule = hp2_balance4)

def hp2_balance5(m, T):
    if not T == m.T.first():
        return m.bSwitchOn_hp2[T] >= m.bOn_hp2[T] - m.bOn_hp2[T-1]
    else:
        return Constraint.Skip
model.hp2_bal_cons5 = Constraint(model.T, rule = hp2_balance5)

def hp2_balance6(m, T):
    if not T == m.T.first():
        return m.bSwitchOn_hp2[T] + m.bSwitchOn_hp2[T - 1] <= 1
    else:
        return Constraint.Skip
model.hp2_bal_cons6 = Constraint(model.T, rule = hp2_balance6)

def hp3_balance1(m, T):
    return m.P_th_cool_hp3[T] + m.P_el_hp3[T] == m.P_th_heat_hp3[T]
model.hp3_bal_cons1 = Constraint(model.T, rule = hp3_balance1)

def hp3_balance2(m, T):
    return ((m.T_flow_heat[T] * m.corr_fac_hp)/(m.T_flow_heat[T] - m.T_flow_cool[T])) * m.P_el_hp3[T] - m.bSwitchOn_hp3[T] * m.switch_loss_hp * 125 == m.P_th_heat_hp3[T]
model.hp3_balance2 = Constraint(model.T, rule = hp3_balance2)

def hp3_balance3(m, T):
    return m.bOn_hp3[T] * 125 >= m.P_el_hp3[T]
model.hp3_bal_cons3 = Constraint(model.T, rule = hp3_balance3)

def hp3_balance4(m, T):
    return m.P_el_hp3[T] >= m.bOn_hp3[T] * m.rel_min_hp * 125
model.hp3_bal_cons4 = Constraint(model.T, rule = hp3_balance4)

def hp3_balance5(m, T):
    if not T == m.T.first():
        return m.bSwitchOn_hp3[T] >= m.bOn_hp3[T] - m.bOn_hp3[T-1]
    else:
        return Constraint.Skip
model.hp3_bal_cons5 = Constraint(model.T, rule = hp3_balance5)

def hp3_balance6(m, T):
    if not T == m.T.first():
        return m.bSwitchOn_hp3[T] + m.bSwitchOn_hp3[T - 1] <= 1
    else:
        return Constraint.Skip
model.hp3_bal_cons6 = Constraint(model.T, rule = hp3_balance6)

def chp1_balance1(m, T):
    return m.P_gas_chp1[T] * m.eta_th_chp - m.bSwitchOn_chp1[T] * m.switch_loss_chp * 6000 == m.P_th_heat_chp1[T]
model.chp1_balance1 = Constraint(model.T, rule = chp1_balance1)

def chp1_balance2(m, T):
    return m.P_gas_chp1[T] * m.eta_el_chp - m.bSwitchOn_chp1[T] * m.switch_loss_chp * 6000 == m.P_el_chp1[T]
model.chp1_balance2 = Constraint(model.T, rule = chp1_balance2)

def chp1_balance3(m, T):
    return m.bOn_chp1[T] * 6000 >= m.P_gas_chp1[T]
model.chp1_balance3 = Constraint(model.T, rule = chp1_balance3)

def chp1_balance4(m, T):
    return m.P_gas_chp1[T] >= m.bOn_chp1[T] * m.rel_min_chp * 6000
model.chp1_balance4 = Constraint(model.T, rule = chp1_balance4)

def chp1_balance5(m, T):
    if not T == m.T.first():
        return m.bSwitchOn_chp1[T] >= m.bOn_chp1[T] - m.bOn_chp1[T-1]
    else:
        return Constraint.Skip
model.chp1_balance5 = Constraint(model.T, rule = chp1_balance5)

def chp1_balance6(m, T):
    if not T == m.T.first():
        return m.bSwitchOn_chp1[T] + m.bSwitchOn_chp1[T - 1] <= 1
    else:
        return Constraint.Skip
model.chp1_balance6 = Constraint(model.T, rule = chp1_balance6)

# def chp2_balance1(m, T):
#     return m.P_gas_chp2[T] * m.eta_th_chp - m.bSwitchOn_chp2[T] * m.switch_loss_chp * 3000 == m.P_th_heat_chp2[T]
# model.chp2_balance1 = Constraint(model.T, rule = chp2_balance1)

# def chp2_balance2(m, T):
#     return m.P_gas_chp2[T] * m.eta_el_chp - m.bSwitchOn_chp2[T] * m.switch_loss_chp * 3000 == m.P_el_chp2[T]
# model.chp2_balance2 = Constraint(model.T, rule = chp2_balance2)

# def chp2_balance3(m, T):
#     return m.bOn_chp2[T] * 3000 >= m.P_gas_chp2[T]
# model.chp2_balance3 = Constraint(model.T, rule = chp2_balance3)

# def chp2_balance4(m, T):
#     return m.P_gas_chp2[T] >= m.bOn_chp2[T] * m.rel_min_chp * 3000
# model.chp2_balance4 = Constraint(model.T, rule = chp2_balance4)

# def chp2_balance5(m, T):
#     if not T == m.T.first():
#         return m.bSwitchOn_chp2[T] >= m.bOn_chp2[T] - m.bOn_chp2[T-1]
#     else:
#         return Constraint.Skip
# model.chp2_balance5 = Constraint(model.T, rule = chp2_balance5)

# def chp2_balance6(m, T):
#     if not T == m.T.first():
#         return m.bSwitchOn_chp2[T] + m.bSwitchOn_chp2[T - 1] <= 1
#     else:
#         return Constraint.Skip
# model.chp2_balance6 = Constraint(model.T, rule = chp2_balance6)

def ct_balance1(m, T):
    return (50 - (45/40)*(m.T_amb[T] - 273.15)) * m.P_el_ct[T] == m.P_th_cool_ct[T]
model.ct_balance1 = Constraint(model.T, rule = ct_balance1)

# def influence_factor_def1(m, T):
#     return m.controlLogic.fac[1, T] == m.c_el[T]
# model.infl_fac_def1 = Constraint(model.T, rule = influence_factor_def1)

# def influence_factor_def2(m, T):
#     return m.controlLogic.fac[2, T] == m.T_amb[T]
# model.infl_fac_def2 = Constraint(model.T, rule = influence_factor_def2)

# def influence_factor_def3(m, T):
#     return m.controlLogic.fac[3, T] == m.T_amb[T]
# model.infl_fac_def3 = Constraint(model.T, rule = influence_factor_def3)

# def influence_factor_def4(m, T):
#     return m.controlLogic.fac[4, T] == m.T_amb[T]
# model.infl_fac_def4 = Constraint(model.T, rule = influence_factor_def4)

# def rel_load_def1(m, T):
#     return m.controlLogic.rel[11, T] == m.P_el_hp1[T] / 375
# model.rel_load_def1 = Constraint(model.T, rule = rel_load_def1)

# def rel_load_def2(m, T):
#     return m.controlLogic.rel[22, T] == m.P_el_hp2[T] / 250
# model.rel_load_def2 = Constraint(model.T, rule = rel_load_def2)

# def rel_load_def3(m, T):
#     return m.controlLogic.rel[33, T] == m.P_el_hp3[T] / 125
# model.rel_load_def3 = Constraint(model.T, rule = rel_load_def3)

# def rel_load_def4(m, T):
#     return m.controlLogic.rel[44, T] == m.P_gas_chp1[T] / 6000
# model.rel_load_def4 = Constraint(model.T, rule = rel_load_def4)

# def rel_load_def5(m, T):
#     return m.controlLogic.rel[55, T] == m.P_gas_chp2[T] / 3000
# model.rel_load_def5 = Constraint(model.T, rule = rel_load_def5

def obj(m):
    return sum((m.P_el_hp1[t] + m.P_el_hp2[t] + m.P_el_hp3[t] + m.P_el_ct[t] - m.P_el_chp1[t])*m.c_el[t] + (m.P_gas_chp1[t])*m.c_gas for t in m.T)
model.obj_rule = Objective(rule = obj)

print('data loaded')
opt = SolverFactory('cplex')
instance = model.create_instance(model_data)
print('model set up')
opt.solve(instance, tee = True, options_string="mipgap=0.01")

var_list = [
    'P_th_heat_hp1',
    'P_th_heat_hp2',
    'P_th_heat_hp3',
    'P_th_cool_hp1',
    'P_th_cool_hp2',
    'P_th_cool_hp3',
    'P_el_hp1',
    'P_el_hp2',
    'P_el_hp3',
    'P_th_heat_chp1',
    'P_el_chp1',
    'P_gas_chp1',
    'P_th_cool_ct',
    'P_el_ct'
    ]

param_list = [
    'T_amb',
    'T_amb_avg',
    'T_flow_heat',
    'T_flow_cool',
    'P_th_heat',
    'P_th_cool',
    'c_el'
    ]

res_df = pd.DataFrame()
for var in var_list:
    temp_list = []
    for key in getattr(instance, var)._data:
        temp_list.append(getattr(instance, var)._data[key].value)
        
    res_df[var] = temp_list

for param in param_list:
    temp_list = []
    for key in getattr(instance, param)._data:
        temp_list.append(getattr(instance, param)._data[key])
        
    res_df[param] = temp_list


# for conv in [11, 22, 33, 44]:
#     temp_list = []
#     for timestep in getattr(instance, 'T')._ordered_values:
#         for key in getattr(instance.controlLogic, 'c1')._data:
#             if getattr(instance.controlLogic, 'c1')._data[key].value > 0.8 and key[0] == conv and key[2] == timestep:
#                 temp_list.append(key[1])
#     res_df[str(conv)] = temp_list

res_df.to_excel(os.path.join(current_dir, "res_" + scen + "_" + type + ".xlsx"))








