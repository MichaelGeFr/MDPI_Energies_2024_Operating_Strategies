#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   utilities.py
@Time    :   2023/02/20 13:04:00
@Author  :   Lukas Theisinger 
@Version :   1.0
@Contact :   l.theisinger@ptw.tu-darmstadt.de
'''

import pandas as pd

def prepare_timeseries(path, indexing_dict):
    input_dict = {}
    input_dict[None] = {}

    raw_df = pd.read_excel(path)

    for key in indexing_dict:
        if not indexing_dict[key] == None:
            input_dict[None][key] = {}
            for index, row in raw_df.iterrows():
                input_dict[None][key][row[indexing_dict[key]]] = row[key]
        else:
            input_dict[None][key] = {}
            input_dict[None][key][None] = list(raw_df[key].values)

    return input_dict

def controlLogic_data_c1_optimization(number_of_converters, T):
    input_dict = {}
    input_dict[None] = {}
    input_dict[None]['T'] = {}
    input_dict[None]['T'][None] = T
    input_dict[None]['P'] = {}
    input_dict[None]['P'][None] = [i for i in range(1, number_of_converters + 1)]
    input_dict[None]['C'] = {}
    input_dict[None]['C'][None] = [i*11 for i in range(1, number_of_converters + 1)]
    
    return input_dict

def controlLogic_validation(thresholds, rules, number_of_influence_factors, number_of_converters, T):        

    input_dict = {}
    input_dict[None] = {}
    input_dict[None]['T'] = {}
    input_dict[None]['T'][None] = T
    input_dict[None]['I'] = {}
    input_dict[None]['I'][None] = [i for i in range(1, number_of_influence_factors + 1)]
    input_dict[None]['J'] = {}
    input_dict[None]['J'][None] = [i for i in range(1, pow(2, number_of_influence_factors) + 1)]
    input_dict[None]['P'] = {}
    input_dict[None]['P'][None] = [i for i in range(1, number_of_converters + 1)]
    input_dict[None]['C'] = {}
    input_dict[None]['C'][None] = [i*11 for i in range(1, number_of_converters + 1)]
    input_dict[None]['A'] = {}
    for i in input_dict[None]['I'][None]:
        input_dict[None]['A'][i] = pow(2, i-1)
    input_dict[None]['B'] = {}
    for j in input_dict[None]['J'][None]:
        input_dict[None]['B'][j] = j

    input_dict[None]['c2'] = {}

    for conv in input_dict[None]['C'][None]:
        for prio in input_dict[None]['P'][None]:
            for rule in input_dict[None]['J'][None]:
                if rules[rule][prio - 1] == conv:
                    input_dict[None]['c2'][(conv, prio, rule)] = 1
                else:
                    input_dict[None]['c2'][(conv, prio, rule)] = 0

    input_dict[None]['thres'] = {}

    for i in input_dict[None]['I'][None]:
        input_dict[None]['thres'][i] = thresholds[i]

    return input_dict


