#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   ControlLogic.py
@Time    :   2023/02/17 11:00:00
@Author  :   Lukas Theisinger 
@Version :   1.0
@Contact :   l.theisinger@ptw.tu-darmstadt.de
'''

from pyomo.environ import *

ControlLogic = AbstractModel()

### indices ###
ControlLogic.T = Set(domain = NonNegativeIntegers, doc = "time step index T")
ControlLogic.P = Set(domain = NonNegativeIntegers, doc = "priority index P")
ControlLogic.C = Set(doc = "converter index C")

### variables ###
ControlLogic.c1 = Var(ControlLogic.C, ControlLogic.P, ControlLogic.T, domain = Binary, doc = "fullfillment variable to determine, wheter converter C is at priority P (C, P, T)")
ControlLogic.d = Var(ControlLogic.C, ControlLogic.P, ControlLogic.T, domain = Binary, doc = "approval variable of converter C (C, P, T)")
ControlLogic.fOperatingPoint = Var(ControlLogic.C, ControlLogic.P, ControlLogic.T, bounds = (0, 1), doc = "relative utilization of converter C (C, P, T)")
ControlLogic.rel = Var(ControlLogic.C, ControlLogic.T, bounds = (0, 1), doc = "relative utilization of converter C (C, T)")

### parameters ###
ControlLogic.M = 100


### constraints ###
def single_converter_constraint(m, P, T):
    """
    per priority only one converter 
    """
    return sum(m.c1[c, P, T] for c in m.C) == 1
ControlLogic.single_conv_cons = Constraint(ControlLogic.P, ControlLogic.T, rule = single_converter_constraint)

def single_priority_constraint(m, C, T):
    """
    converter only at one priority 
    """
    return sum(m.c1[C, p, T] for p in m.P) == 1
ControlLogic.single_prio_cons = Constraint(ControlLogic.C, ControlLogic.T, rule = single_priority_constraint)

def approval_constraint1(m, C, P, T):
    """ 
    if fOperatingPoint of prior converter is 1 then d = 1 else d = 0
    """
    if P == m.P.first():
        return sum(m.d[c, P, T] for c in m.C) == 1
    else:
        return Constraint.Skip
ControlLogic.approval_cons1 = Constraint(ControlLogic.C, ControlLogic.P, ControlLogic.T, rule = approval_constraint1)

def approval_constraint2(m, C, P, T):
    """ 
    if fOperatingPoint of prior converter is 1 then d = 1 else d = 0
    """
    if not P == m.P.first():
        return m.d[C, P, T] <= sum(m.fOperatingPoint[c, P-1, T] for c in m.C)
    else:
        return Constraint.Skip
ControlLogic.approval_cons2 = Constraint(ControlLogic.C, ControlLogic.P, ControlLogic.T, rule = approval_constraint2)

def approval_constraint3(m, C, P, T):
    """ 
    allow approval only if priority is present
    """
    return m.d[C, P, T] <= m.c1[C, P, T]
ControlLogic.approval_cons3 = Constraint(ControlLogic.C, ControlLogic.P, ControlLogic.T, rule = approval_constraint3)

def part_load_constraint1(m, C, P, T):
    """
    definition operating point
    """
    return sum(m.fOperatingPoint[C, p, T] for p in m.P) == m.rel[C, T]
ControlLogic.partl_cons1 = Constraint(ControlLogic.C, ControlLogic.P, ControlLogic.T, rule = part_load_constraint1)

def part_load_constraint2(m, C, P, T):
    """
    enforce off-status without activation
    """
    return m.fOperatingPoint[C, P, T] <= m.d[C, P, T]
ControlLogic.partl_cons2 = Constraint(ControlLogic.C, ControlLogic.P, ControlLogic.T, rule = part_load_constraint2)