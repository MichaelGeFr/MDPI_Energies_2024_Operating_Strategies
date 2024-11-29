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
ControlLogic.I = Set(domain = NonNegativeIntegers, doc = "influence factor index I")
ControlLogic.J = Set(domain = NonNegativeIntegers, doc = "decision rule index J")
ControlLogic.T = Set(domain = NonNegativeIntegers, doc = "time step index T")
ControlLogic.P = Set(domain = NonNegativeIntegers, doc = "priority index P")
ControlLogic.C = Set(doc = "converter index C")

### variables ###
ControlLogic.fac = Var(ControlLogic.I, ControlLogic.T, domain = Reals, doc = "influence factor (I, T)")
ControlLogic.thres = Var(ControlLogic.I, domain = Reals, doc = "threshold values (I)")

ControlLogic.a = Var(ControlLogic.I, ControlLogic.T, domain = Binary, doc = "fullfillment variable for combination of influence factor and threshold value (I, T)")
ControlLogic.b = Var(ControlLogic.J, ControlLogic.T, domain = Binary, doc = "decision variable to determine the activation of a decision rule (J, T)")
ControlLogic.c1 = Var(ControlLogic.C, ControlLogic.P, ControlLogic.T, domain = Binary, doc = "fullfillment variable to determine, wheter converter C is at priority P (C, P, T)")
ControlLogic.c2 = Var(ControlLogic.C, ControlLogic.P, ControlLogic.J, domain = Binary, doc = "fullfillment variable to determine, wheter converter C is at priority P in decision rule J (C, P, J)") 
ControlLogic.c3 = Var(ControlLogic.C, ControlLogic.P, ControlLogic.J, ControlLogic.T, domain = Binary, doc = "fullfillment variable to determine, wheter converter C is at priority P in decision rule J at time T (C, P, J ,T)") 
ControlLogic.d = Var(ControlLogic.C, ControlLogic.P, ControlLogic.T, domain = Binary, doc = "approval variable of converter C (C, P, T)")
ControlLogic.less = Var(ControlLogic.I, ControlLogic.T, domain = Binary, doc = "factor less than threshold (I, T)")
ControlLogic.greather = Var(ControlLogic.I, ControlLogic.T, domain = Binary, doc = "factor greather than threshold (I, T)")

ControlLogic.fOperatingPoint = Var(ControlLogic.C, ControlLogic.P, ControlLogic.T, domain = NonNegativeReals, doc = "relative utilization of converter C (C, P, T)")
ControlLogic.bSetStatusOn = Var(ControlLogic.C, ControlLogic.P, ControlLogic.T, domain = Binary, doc = "activation of converter C (C, P, T)")
ControlLogic.rel = Var(ControlLogic.C, ControlLogic.T, domain = NonNegativeReals, doc = "relative utilization of converter C (C, T)")

### parameters ###
ControlLogic.A = Param(ControlLogic.I, doc = "[1, 2, 4.. 2^i]")
ControlLogic.B = Param(ControlLogic.J, doc = "[1, 2, 3.. j]")
ControlLogic.M = 100

ControlLogic.P_nom = Param(ControlLogic.C)


### constraints ###

def influence_threshold_constraint1(m, I, T):
    """
    if fac >= thres then a = 1 else a = 0
    """
    return m.fac[I, T] + (1 - m.a[I, T]) * m.M >= m.thres[I]
ControlLogic.infl_thres_cons1 = Constraint(ControlLogic.I, ControlLogic.T, rule = influence_threshold_constraint1)

def influence_threshold_constraint2(m, I, T):
    """
    if fac >= thres then a = 1 else a = 0
    """
    return m.fac[I, T] - m.a[I, T] * m.M <= m.thres[I]
ControlLogic.infl_thres_cons2 = Constraint(ControlLogic.I, ControlLogic.T, rule = influence_threshold_constraint2)

def influence_threshold_constraint3(m, I, T):
    """
    do not allow equality of fac and thres
    """
    return m.fac[I, T] + m.less[I, T] * m.M >= m.thres[I] + (0.1) # adjusted
ControlLogic.infl_thres_cons3 = Constraint(ControlLogic.I, ControlLogic.T, rule = influence_threshold_constraint3)

def influence_threshold_inequality1(m, I, T):
    """
    do not allow equality of fac and thres
    """
    return m.fac[I, T] - m.greather[I, T] * m.M <= m.thres[I] - (0.1) # adjusted
ControlLogic.infl_thres_ineq1 = Constraint(ControlLogic.I, ControlLogic.T, rule = influence_threshold_inequality1)


def influence_threshold_inequality2(m, I, T):
    """
    do not allow equality of fac and thres
    """
    return m.greather[I, T] + m.less[I, T] == 1
ControlLogic.infl_thres_ineq2 = Constraint(ControlLogic.I, ControlLogic.T, rule = influence_threshold_inequality2)


def decision_rule_constraint1(m, J, T):
    """
    linkage between threshold-activation and decision rules
    """
    return sum(m.b[j,T] * m.B[j] for j in m.J) == sum(m.a[i, T] * m.A[i] for i in m.I) + 1
ControlLogic.dec_rule_cons1 = Constraint(ControlLogic.J, ControlLogic.T, rule = decision_rule_constraint1)

def decision_rule_constraint2(m,T):
    """
    only one decision rule active
    """
    return sum(m.b[j, T] for j in m.J) == 1
ControlLogic.dec_rule_cons2 = Constraint(ControlLogic.T, rule = decision_rule_constraint2)

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

def single_converter_constraint2(m, P, J):
    """
    per priority only one converter
    """
    return sum(m.c2[c, P, J] for c in m.C) == 1
ControlLogic.single_conv_cons2 = Constraint(ControlLogic.P, ControlLogic.J, rule = single_converter_constraint2)

def single_priority_constraint2(m, C, J):
    """
    converter only at one priority
    """
    return sum(m.c2[C, p, J] for p in m.P) == 1
ControlLogic.single_prio_cons2 = Constraint(ControlLogic.C, ControlLogic.J, rule = single_priority_constraint2)

def prio_selection_constraint1(m, C, P, J, T):
    """ 
    if decision rule not active then c3 = 0
    """
    return m.b[J, T] >= m.c3[C, P, J, T]
ControlLogic.prio_selec_cons1 = Constraint(ControlLogic.C, ControlLogic.P, ControlLogic.J, ControlLogic.T, rule = prio_selection_constraint1)

def prio_selection_constraint2(m, C, P, J, T):
    """ 
    if decision rule active then c3 <= c2
    """
    return m.c3[C, P, J, T] <= m.c2[C, P, J] + (1 - m.b[J, T]) * m.M
ControlLogic.prio_selec_cons2 = Constraint(ControlLogic.C, ControlLogic.P, ControlLogic.J, ControlLogic.T, rule = prio_selection_constraint2)

def prio_selection_constraint3(m, C, P, J, T):
    """ 
    if decision rule active then c3 >= c2
    """
    return (1 - m.b[J, T]) * m.M + m.c3[C, P, J, T] >= m.c2[C, P, J]
ControlLogic.prio_selec_cons3 = Constraint(ControlLogic.C, ControlLogic.P, ControlLogic.J, ControlLogic.T, rule = prio_selection_constraint3)

def prio_selection_constraint4(m, C, P, T):
    """ 
    enforce priority of decision rules
    """
    return m.c1[C, P, T] == sum(m.c3[C, P, j, T] for j in m.J)
ControlLogic.prio_selec_cons4 = Constraint(ControlLogic.C, ControlLogic.P, ControlLogic.T, rule = prio_selection_constraint4)

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

def activation_constraint1(m, C, P, T):
    """
    allow activation only if approval is present
    """
    return m.d[C, P, T] >= m.bSetStatusOn[C, P, T]
ControlLogic.act_cons = Constraint(ControlLogic.C, ControlLogic.P, ControlLogic.T, rule = activation_constraint1)

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
    return m.fOperatingPoint[C, P, T] <= m.bSetStatusOn[C, P, T]
ControlLogic.partl_cons2 = Constraint(ControlLogic.C, ControlLogic.P, ControlLogic.T, rule = part_load_constraint2)