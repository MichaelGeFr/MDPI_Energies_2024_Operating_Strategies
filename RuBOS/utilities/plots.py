import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import numpy as np



###
winter_flex = pd.read_excel("C://Users//L.Theisinger_lokal//Documents//GitHub//PTW_git//dissertation_lt//res_winter_flex.xlsx")
winter_sequencing = pd.read_excel("C://Users//L.Theisinger_lokal//Documents//GitHub//PTW_git//dissertation_lt//res_winter_sequencing.xlsx")
winter_rb = pd.read_excel("C://Users//L.Theisinger_lokal//Documents//GitHub//PTW_git//dissertation_lt//res_winter_rb.xlsx")

spring_flex = pd.read_excel("C://Users//L.Theisinger_lokal//Documents//GitHub//PTW_git//dissertation_lt//res_spring_flex.xlsx")
spring_sequencing = pd.read_excel("C://Users//L.Theisinger_lokal//Documents//GitHub//PTW_git//dissertation_lt//res_spring_sequencing.xlsx")
spring_rb = pd.read_excel("C://Users//L.Theisinger_lokal//Documents//GitHub//PTW_git//dissertation_lt//res_spring_rb.xlsx")

summer_flex = pd.read_excel("C://Users//L.Theisinger_lokal//Documents//GitHub//PTW_git//dissertation_lt//res_summer_flex.xlsx")
summer_sequencing = pd.read_excel("C://Users//L.Theisinger_lokal//Documents//GitHub//PTW_git//dissertation_lt//res_summer_sequencing.xlsx")
summer_rb = pd.read_excel("C://Users//L.Theisinger_lokal//Documents//GitHub//PTW_git//dissertation_lt//res_summer_rb.xlsx")


norm = Normalize(vmin=0, vmax=1)
cmap = cm.Greys
m = cm.ScalarMappable(norm=norm, cmap=cmap)

winter_chp1_flex = [m.to_rgba(x) for x in winter_flex['P_gas_chp1']/6000]
winter_hp1_flex = [m.to_rgba(x) for x in winter_flex['P_el_hp1']/375]
winter_hp2_flex = [m.to_rgba(x) for x in winter_flex['P_el_hp2']/250]
winter_hp3_flex = [m.to_rgba(x) for x in winter_flex['P_el_hp3']/125]

winter_chp1_sequencing = [m.to_rgba(x) for x in winter_sequencing['P_gas_chp1']/6000]
winter_hp1_sequencing = [m.to_rgba(x) for x in winter_sequencing['P_el_hp1']/375]
winter_hp2_sequencing = [m.to_rgba(x) for x in winter_sequencing['P_el_hp2']/250]
winter_hp3_sequencing = [m.to_rgba(x) for x in winter_sequencing['P_el_hp3']/125]

winter_chp1_rb = [m.to_rgba(x) for x in winter_rb['P_gas_chp1']/6000]
winter_hp1_rb = [m.to_rgba(x) for x in winter_rb['P_el_hp1']/375]
winter_hp2_rb = [m.to_rgba(x) for x in winter_rb['P_el_hp2']/250]
winter_hp3_rb = [m.to_rgba(x) for x in winter_rb['P_el_hp3']/125]

spring_chp1_flex = [m.to_rgba(x) for x in spring_flex['P_gas_chp1']/6000]
spring_hp1_flex = [m.to_rgba(x) for x in spring_flex['P_el_hp1']/375]
spring_hp2_flex = [m.to_rgba(x) for x in spring_flex['P_el_hp2']/250]
spring_hp3_flex = [m.to_rgba(x) for x in spring_flex['P_el_hp3']/125]

spring_chp1_sequencing = [m.to_rgba(x) for x in spring_sequencing['P_gas_chp1']/6000]
spring_hp1_sequencing = [m.to_rgba(x) for x in spring_sequencing['P_el_hp1']/375]
spring_hp2_sequencing = [m.to_rgba(x) for x in spring_sequencing['P_el_hp2']/250]
spring_hp3_sequencing = [m.to_rgba(x) for x in spring_sequencing['P_el_hp3']/125]

spring_chp1_rb = [m.to_rgba(x) for x in spring_rb['P_gas_chp1']/6000]
spring_hp1_rb = [m.to_rgba(x) for x in spring_rb['P_el_hp1']/375]
spring_hp2_rb = [m.to_rgba(x) for x in spring_rb['P_el_hp2']/250]
spring_hp3_rb = [m.to_rgba(x) for x in spring_rb['P_el_hp3']/125]

summer_chp1_flex = [m.to_rgba(x) for x in summer_flex['P_gas_chp1']/6000]
summer_hp1_flex = [m.to_rgba(x) for x in summer_flex['P_el_hp1']/375]
summer_hp2_flex = [m.to_rgba(x) for x in summer_flex['P_el_hp2']/250]
summer_hp3_flex = [m.to_rgba(x) for x in summer_flex['P_el_hp3']/125]

summer_chp1_sequencing = [m.to_rgba(x) for x in summer_sequencing['P_gas_chp1']/6000]
summer_hp1_sequencing = [m.to_rgba(x) for x in summer_sequencing['P_el_hp1']/375]
summer_hp2_sequencing = [m.to_rgba(x) for x in summer_sequencing['P_el_hp2']/250]
summer_hp3_sequencing = [m.to_rgba(x) for x in summer_sequencing['P_el_hp3']/125]

summer_chp1_rb = [m.to_rgba(x) for x in summer_rb['P_gas_chp1']/6000]
summer_hp1_rb = [m.to_rgba(x) for x in summer_rb['P_el_hp1']/375]
summer_hp2_rb = [m.to_rgba(x) for x in summer_rb['P_el_hp2']/250]
summer_hp3_rb = [m.to_rgba(x) for x in summer_rb['P_el_hp3']/125]

### 
plt.rcParams["font.family"] = "sans-serif"
#plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rcParams.update({'font.size': 10})
plt.rcParams['figure.dpi'] = 100

px = 1/plt.rcParams['figure.dpi']
cm = 1/2.54
fig, ax = plt.subplots(5, 3, figsize=(17*cm, 20*cm))
plt.subplots_adjust(wspace=0.02, hspace=0.4)

ax[0,0].title.set_text('Winter')
ax[0,0].set_ylabel('Power in MW', labelpad = 6.0)
l1, = ax[0,0].plot(winter_flex.index, winter_flex['P_th_heat']/1000, color = "black", label = "Heating demand")
ax[0,0].set_yticks([0, 5], [0.0, 5.0], fontsize = 8)
ax[0,0].set_xticks([])
ax_temp = ax[0,0].twinx()
l2, = ax_temp.plot(winter_flex.index, winter_flex['T_amb'] - 273.15, color = "grey", label = "Ambient temperature")
ax_temp.set_ylim([-5, 35])
ax_temp.set_yticks([])
ax_temp.set_xticks([])

ax[1,0].set_ylabel('MPC', labelpad = 2.0)
ax[1,0].bar(winter_flex.index, height = 1.0, width = 1.0, bottom = 6.0, color = winter_chp1_flex)
ax[1,0].bar(winter_flex.index, height = 1.0, width = 1.0, bottom = 4.0, color = winter_hp1_flex)
ax[1,0].bar(winter_flex.index, height = 1.0, width = 1.0, bottom = 2.0, color = winter_hp2_flex)
ax[1,0].bar(winter_flex.index, height = 1.0, width = 1.0, bottom = 0.0, color = winter_hp3_flex)
ax[1,0].set_xticks([])
ax[1,0].set_yticks(ticks=[0.5, 2.5, 4.5, 6.5], labels=["HP3", "HP2", "HP1", "CHP"], fontsize = 8)

ax[2,0].set_ylabel('Sequencing', labelpad = 2.0)
ax[2,0].bar(winter_flex.index, height = 1.0, width = 1.0, bottom = 6.0, color = winter_chp1_sequencing)
ax[2,0].bar(winter_flex.index, height = 1.0, width = 1.0, bottom = 4.0, color = winter_hp1_sequencing)
ax[2,0].bar(winter_flex.index, height = 1.0, width = 1.0, bottom = 2.0, color = winter_hp2_sequencing)
ax[2,0].bar(winter_flex.index, height = 1.0, width = 1.0, bottom = 0.0, color = winter_hp3_sequencing)
ax[2,0].set_xticks([])
ax[2,0].set_yticks(ticks=[0.5, 2.5, 4.5, 6.5], labels=["HP3", "HP2", "HP1", "CHP"], fontsize = 8)

ax[3,0].set_ylabel('Baseline', labelpad = 2.0)
ax[3,0].set_xlabel('Time in d')
ax[3,0].bar(winter_flex.index, height = 1.0, width = 1.0, bottom = 6.0, color = winter_chp1_rb)
ax[3,0].bar(winter_flex.index, height = 1.0, width = 1.0, bottom = 4.0, color = winter_hp1_rb)
ax[3,0].bar(winter_flex.index, height = 1.0, width = 1.0, bottom = 2.0, color = winter_hp2_rb)
ax[3,0].bar(winter_flex.index, height = 1.0, width = 1.0, bottom = 0.0, color = winter_hp3_rb)
ax[3,0].set_yticks(ticks=[0.5, 2.5, 4.5, 6.5], labels=["HP3", "HP2", "HP1", "CHP"], fontsize = 8)
ax[3,0].set_xticks(ticks=[24, 72, 120, 168], labels=[1, 3, 5, 7])

# cmap = cm.Greys
# norm = mpl.colors.Normalize(vmin=0, vmax=1)
cb1 = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax[4,1], orientation='horizontal', aspect = 16)
cb1.set_ticks([0,1])
cb1.set_ticklabels(["off", "full-load"])


ax[4,1].legend(handles = [l1, l2], loc="upper center", bbox_to_anchor=(0.5, 0.75), frameon = False)
ax[4,0].axis('off')
ax[4,1].axis('off')
ax[4,2].axis('off')
# ax[5,0].axis('off')
# ax[5,1].axis('off')
# ax[5,2].axis('off')


ax[0,1].title.set_text('Spring')
ax[0,1].plot(spring_flex.index, spring_flex['P_th_heat']/1000, color = "black")
ax[0,1].set_ylim([0, 5])
ax[0,1].set_xticks([])
ax[0,1].set_yticks([])
ax_temp = ax[0,1].twinx()
ax_temp.plot(spring_flex.index, spring_flex['T_amb'] - 273.15, color = "grey")
ax_temp.set_ylim([-5, 35])
ax_temp.set_xticks([])
ax_temp.set_yticks([])

ax[1,1].bar(spring_flex.index, height = 1.0, width = 1.0, bottom = 6.0, color = spring_chp1_flex)
ax[1,1].bar(spring_flex.index, height = 1.0, width = 1.0, bottom = 4.0, color = spring_hp1_flex)
ax[1,1].bar(spring_flex.index, height = 1.0, width = 1.0, bottom = 2.0, color = spring_hp2_flex)
ax[1,1].bar(spring_flex.index, height = 1.0, width = 1.0, bottom = 0.0, color = spring_hp3_flex)
ax[1,1].set_xticks([])
ax[1,1].set_yticks([])


ax[2,1].bar(spring_flex.index, height = 1.0, width = 1.0, bottom = 6.0, color = spring_chp1_sequencing)
ax[2,1].bar(spring_flex.index, height = 1.0, width = 1.0, bottom = 4.0, color = spring_hp1_sequencing)
ax[2,1].bar(spring_flex.index, height = 1.0, width = 1.0, bottom = 2.0, color = spring_hp2_sequencing)
ax[2,1].bar(spring_flex.index, height = 1.0, width = 1.0, bottom = 0.0, color = spring_hp3_sequencing)
ax[2,1].set_xticks([])
ax[2,1].set_yticks([])

ax[3,1].set_xlabel('Time in d')
ax[3,1].bar(spring_flex.index, height = 1.0, width = 1.0, bottom = 6.0, color = spring_chp1_rb)
ax[3,1].bar(spring_flex.index, height = 1.0, width = 1.0, bottom = 4.0, color = spring_hp1_rb)
ax[3,1].bar(spring_flex.index, height = 1.0, width = 1.0, bottom = 2.0, color = spring_hp2_rb)
ax[3,1].bar(spring_flex.index, height = 1.0, width = 1.0, bottom = 0.0, color = spring_hp3_rb)
ax[3,1].set_yticks([])
ax[3,1].set_xticks(ticks=[24, 72, 120, 168], labels=[1, 3, 5, 7])

ax[0,2].title.set_text('Summer')
ax[0,2].plot(summer_flex.index, summer_flex['P_th_heat']/1000, color = "black")
ax[0,2].set_ylim([0, 5])
ax[0,2].set_xticks([])
ax[0,2].set_yticks([])
ax_temp = ax[0,2].twinx()
ax_temp.plot(summer_flex.index, summer_flex['T_amb'] - 273.15, color = "grey")
ax_temp.set_ylabel('Temperature in °C', labelpad = 1.0)
ax_temp.set_yticks([-5, 35], [-5, 35], fontsize = 8)
ax_temp.set_ylim([-5, 35])
ax_temp.set_xticks([])

ax[1,2].bar(summer_flex.index, height = 1.0, width = 1.0, bottom = 6.0, color = summer_chp1_flex)
ax[1,2].bar(summer_flex.index, height = 1.0, width = 1.0, bottom = 4.0, color = summer_hp1_flex)
ax[1,2].bar(summer_flex.index, height = 1.0, width = 1.0, bottom = 2.0, color = summer_hp2_flex)
ax[1,2].bar(summer_flex.index, height = 1.0, width = 1.0, bottom = 0.0, color = summer_hp3_flex)
ax[1,2].set_xticks([])
ax[1,2].set_yticks([])

ax[2,2].bar(summer_flex.index, height = 1.0, width = 1.0, bottom = 6.0, color = summer_chp1_sequencing)
ax[2,2].bar(summer_flex.index, height = 1.0, width = 1.0, bottom = 4.0, color = summer_hp1_sequencing)
ax[2,2].bar(summer_flex.index, height = 1.0, width = 1.0, bottom = 2.0, color = summer_hp2_sequencing)
ax[2,2].bar(summer_flex.index, height = 1.0, width = 1.0, bottom = 0.0, color = summer_hp3_sequencing)
ax[2,2].set_xticks([])
ax[2,2].set_yticks([])

ax[3,2].set_xlabel('Time in d')
ax[3,2].bar(summer_flex.index, height = 1.0, width = 1.0, bottom = 6.0, color = summer_chp1_rb)
ax[3,2].bar(summer_flex.index, height = 1.0, width = 1.0, bottom = 4.0, color = summer_hp1_rb)
ax[3,2].bar(summer_flex.index, height = 1.0, width = 1.0, bottom = 2.0, color = summer_hp2_rb)
ax[3,2].bar(summer_flex.index, height = 1.0, width = 1.0, bottom = 0.0, color = summer_hp3_rb)
ax[3,2].set_yticks([])
ax[3,2].set_xticks(ticks=[24, 72, 120, 168], labels=[1, 3, 5, 7])

plt.show()
