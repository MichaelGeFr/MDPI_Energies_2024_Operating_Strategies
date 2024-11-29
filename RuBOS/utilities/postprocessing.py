import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

res_df = pd.read_excel("C://Users//L.Theisinger_lokal//Documents//GitHub//PTW_git//dissertation_lt//EnergyInfo_data_Merck_identification.xlsx")

prio_sets = []
for index, row in res_df.iterrows():
    
    prio_set = [row['_11'], row['_22'], row['_33'], row['_44']]
    prio_sets.append(prio_set)

unique_prio_sets = [list(x) for x in set(tuple(x) for x in prio_sets)]

prio_sets_dict = {}
counter = 1
for set in unique_prio_sets:
    prio_sets_dict[str(set)] = counter
    counter += 1

res_df['decision_rule'] = pd.Series(np.zeros(len(res_df.index._values)), index=res_df.index)



for index, row in res_df.iterrows():
    res_df.at[index, 'decision_rule'] = prio_sets_dict[str([row['_11'], row['_22'], row['_33'], row['_44']])]

# ### 
# plt.rcParams["font.family"] = "sans-serif"
# plt.rcParams.update({'font.size': 10})
# plt.rcParams['figure.dpi'] = 100

# px = 1/plt.rcParams['figure.dpi']
# cm = 1/2.54
# fig = plt.figure(figsize=(8.5*cm, 6*cm))

# plt.hist(res_df['decision_rule'], bins = 15, color="k")
# plt.xlabel("Priority set", fontsize=14)
# plt.ylabel("Incidence", fontsize=14)
# plt.tight_layout(w_pad=0.1, pad = 0.1)
# plt.show()
# keep = [1,3,5,7,10,11]

# for rule in res_df['decision_rule'].unique():
#     if rule not in keep:
#         res_df.drop(res_df.index[res_df['decision_rule'] == rule], inplace = True)

# color_map = {
#     1: 'r',
#     3: 'g',
#     5: 'r',
#     7: 'g',
#     10: 'b',
#     11: 'b'
# }

markers = {
    1: 'o',
    3: '^',
    5: 'o',
    7: '^',
    10: 'P',
    11: 'P'
}

# color_map = {
#     1: 'r',
#     3: 'g',
#     5: 'c',
#     7: 'm',
#     10: 'y',
#     11: 'b'
# }

circle = res_df.loc[(res_df['decision_rule'] == 1) | (res_df['decision_rule'] == 5)]
up = res_df.loc[(res_df['decision_rule'] == 3) | (res_df['decision_rule'] == 7)]
plus = res_df.loc[(res_df['decision_rule'] == 10) | (res_df['decision_rule'] == 11)]

# ### 
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams.update({'font.size': 10})
plt.rcParams['figure.dpi'] = 100

px = 1/plt.rcParams['figure.dpi']
cm = 1/2.54
fig = plt.figure(figsize=(10*cm, 6*cm))

plt.scatter(x=circle['T_amb'] - 273.15, y=circle['c_el'], marker="o", color = "k")
plt.scatter(x=plus['T_amb'] - 273.15, y=plus['c_el'], marker="P", color = "lightgray")
plt.scatter(x=up['T_amb'] - 273.15, y=up['c_el'], marker="^", color = "gray")
plt.xlabel("Ambient temperature in °C", fontsize=14)
plt.ylabel("Price in €/kWh", fontsize=14)
plt.ylim((0.1,0.3))
#plt.tight_layout(w_pad=0.1, pad = 0.1)
plt.tight_layout()
plt.show()